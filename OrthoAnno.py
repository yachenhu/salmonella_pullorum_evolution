#!/usr/bin/python3

"""Annotate mutation of orthologs in query genome by refering a CDS database."""

import argparse
import sys
import csv
import re
from io import StringIO
import tempfile

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.SeqUtils import seq3
import pybedtools

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--query', '-q', required=True, type=str, metavar='FILENAME',
        help="Query genome filename"
        )
    parser.add_argument(
        '--db', '-d', required=True, type=str, metavar='FILENAME',
        help="Reference CDS database filename"
        )
    parser.add_argument(
        '--ident', '-i', type=float, default=0.95, metavar='FLOAT',
        help="Identity threshold for saving hsps [0.95]"
        )
    parser.add_argument(
        '--evalue', '-e', type=float, default=1e-10, metavar='FLOAT',
        help="Expectation value threshold for saving hsps [1e-10]"
        )
    parser.add_argument(
        '--output', '-o', type=argparse.FileType('w'), default=sys.stdout,
        metavar='FILENAME', help="Output filename"
        )
    return parser.parse_args()

def tblastn(query, subject_filename):
    """tblastn wrapper"""
    with tempfile.NamedTemporaryFile(delete=True) as query_file:
        SeqIO.write(query, query_file.name, 'fasta')
        proc = NcbitblastnCommandline(
            query=query_file.name, subject=subject_filename, 
            outfmt=5, parse_deflines=True, seg='no', db_gencode=11
            )
        out, err = proc()
    if err:
        raise Exception(err)
    else:
        return StringIO(out)

def filter_hsps(blast_record, ident=0, evalue=10):
    """Filter HSP objects in Blast Record and return as an iterator."""
    for hit in blast_record.alignments:
        for hsp in hit.hsps:
            if (hsp.identities / hsp.align_length >= ident and 
                    hsp.expect <= evalue):
                # Inherit attributes to HSP.
                #hsp.query_length = blast_record.query_length
                hsp.sbjct_id = hit.hit_id
                hsp.sbjct_frame = hsp.frame[1]
                hsp.sbjct_strand = '+' if hsp.sbjct_frame > 0 else '-'
                yield hsp
            else:
                pass

def match_length(pattern, string):
    """Return the length of string matched, 0 for no matched"""
    search = re.search(pattern, string)
    if search:
        return len(search.group())
    else:
        return 0

def trim_seq(seq, left, right):
    return(seq[left:len(seq)-right])

def trim_hsp(hsp):
    """Trim ends with inconsistent bases in a HSP object.
    Modify only position attributes of
        query_start, query_end, sbjct_start, sbjct_end.
    Preserve the others
    """
    # Number of bases (' ' or '+') to trim from left
    left = match_length('\A[+ ]+', hsp.match)
    # Number of bases (' ' or '+') to trim from right
    right = match_length('[+ ]+\Z', hsp.match)
    hsp.query_start += left
    hsp.query_end -= right
    if hsp.frame[1] > 0:
        # Match the forward strand
        hsp.sbjct_start += left * 3
        hsp.sbjct_end -= right * 3
    else:
        # Match the reverse strand
        hsp.sbjct_start += right * 3
        hsp.sbjct_end -= left * 3
    return hsp

def hsp_to_interval(hsp):
    chrom = hsp.sbjct_id
    start, end = sorted([hsp.sbjct_start, hsp.sbjct_end], key=int)
    start -= 1
    strand = hsp.sbjct_strand
    return pybedtools.Interval(chrom, start, end, strand=strand)

def hsps_to_bed(hsps):
    if hsps:
        bed = pybedtools.BedTool([hsp_to_interval(hsp) for hsp in hsps])
        return bed.sort().merge(d=150, s=True)
    else:
        return pybedtools.BedTool([])

def bed_to_location(bed):
    """Convert merged BedTool to string of locations"""
    if bed:
        locations = []
        for interval in bed:
            location = '{0}|{1}..{2}|{3}'.format(
                interval.chrom, interval.start+1, interval.end, interval[3]
                )
            locations.append(location)
        return ','.join(locations)
    else:
        return '.'

def len_str(s, ex=[]):
    """Count the length of a string without sub-strings or characters."""
    num_ex = 0
    for i in set(ex):
        num_ex += s.count(i)
    return len(s) - num_ex

def main():
    args = parse_args()
    fields = ['#id', 'annotation', 'variant', 'location', 'coverage']
    writer = csv.DictWriter(args.output, delimiter='\t', fieldnames=fields)
    #writer.writeheader()
    for seq_record in SeqIO.parse(args.db, 'fasta'):
        out = tblastn(seq_record, args.query)
        blast_record = NCBIXML.read(out)
        query = str(seq_record.seq)
        # New empty dict of output row
        row = dict.fromkeys(fields, '.')
        row['#id'] = seq_record.id
        hsps = list(filter_hsps(blast_record, args.ident, args.evalue))
        # Trim hsps.
        hsps = [trim_hsp(hsp) for hsp in hsps]
        # Sort hsps.
        hsps = sorted(hsps, key=lambda x: x.query_start)
        # Convert HSPs to BedTool
        bed = hsps_to_bed(hsps)
        # Extract locations from BedTool
        row['location'] = bed_to_location(bed)
        if len(bed) > 1:
            # Interrupted or duplicated. Need manual inspection.
            row['annotation'] = 'interrupted'
        elif len(bed) == 0:
            # Deleted.
            row['annotation'] = 'deleted'
            row['variant'] = 'p.0?'
            row['coverage'] = 0
        else:
            # Pop first hsp.
            hsp = hsps.pop(0)
            stop = hsp.sbjct.find('*')
            if hsp.query_start != 1:
                # Start lost.
                row['annotation'] = 'start_lost'
                row['variant'] = 'p.{0}1?'.format(seq3(query[0]))
                row['coverage'] = 0
            elif stop > 0 and stop < len(hsp.sbjct)-1:
                # Prematur termination codon.
                row['annotation'] = 'PTC'
                row['variant'] = 'p.{0}{1}*'.format(
                    seq3(query[stop]), stop+1
                    )
                row['coverage'] = '{0:.2f}'.format(
                    len_str(hsp.match[:stop], ['+', ' ', '*'])
                    / len_str(query, ['*']) * 100
                    )
            elif len(query) > hsp.query_end:
                # Frameshift or truncated.
                row['annotation'] = 'frameshift'
                row['variant'] = 'p.{0}{1}fs'.format(
                    seq3(query[hsp.query_end]), hsp.query_end+1
                    )
                row['coverage'] = '{0:.2f}'.format(
                    len_str(hsp.match, ['+', ' ', '*'])
                    / len_str(query, ['*']) * 100
                    )
            else:
                # A intact CDS.
                row['annotation'] = 'intact'
                row['variant'] = 'p.='
                row['coverage'] = '{0:.2f}'.format(
                    len_str(hsp.match, ['+', ' ', '*'])
                    / len_str(query, ['*']) * 100
                    )
        writer.writerow(row)
    args.output.close()

if __name__ == '__main__':
    sys.exit(main())
