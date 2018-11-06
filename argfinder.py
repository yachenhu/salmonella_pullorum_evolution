#!/usr/bin/python3

# Author: Yachen Hu (yachenhu@outlook.com)
# Dependencies:
#     Biopython  http://biopython.org/
#     pybedtools	https://daler.github.io/pybedtools/

"""
Detect acquired antimicrobial resistance genes
"""

import argparse
import os
import sys
import glob
import io
from string import ascii_letters, digits
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML
import csv
import pybedtools

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene_db', required=True, type=argparse.FileType('r'),
                        help="FASTA file of reference orfs")
    parser.add_argument('--input', required=True, type=argparse.FileType('r'),
                        help="FASTA file of genome sequence")
    parser.add_argument('--output', type=argparse.FileType('w'),
                        default=sys.stdout, help="Output file name")
    return parser.parse_args()

def generate_string(size=8, chars=ascii_letters+digits):
    return ''.join(random.sample(chars, k=size))

def clean(prefix):
    for i in glob.glob('{0}.*'.format(prefix)):
        os.remove(i)

def locate(genome, gene_db):
    prefix = generate_string()
    query_fn = '{}.query.fasta'.format(prefix)
    sbjct_fn = '{}.sbjct.fasta'.format(prefix)
    SeqIO.write(genome, query_fn, 'fasta')
    SeqIO.write(gene_db, sbjct_fn, 'fasta')
    blastn_cline = NcbiblastnCommandline(
        query=query_fn, subject=sbjct_fn, evalue=1e-10,
        outfmt="'6 qseqid qstart qend length'"
    )
    stdout, stderr = blastn_cline()
    csv_reader = csv.reader(stdout.splitlines(), delimiter='\t')
    intervals = []
    for row in csv_reader:
        if int(row[-1]) > 100:
            intervals.append([row[0], int(row[1])-1, int(row[2])])
        else:
            pass
    clean(prefix)
    if intervals:
        return pybedtools.BedTool(intervals).sort().merge(d=-100) # deal with overlaps no more than 100 bp
    else:
        return None
            
def get_seq(genome, intervals):
    seq_records = []
    for seq_id, start, end in intervals:
        start = int(start)
        end = int(end)
        feature = SeqFeature(FeatureLocation(start, end))
        seq = feature.extract(SeqIO.to_dict(genome)[seq_id].seq)
        seq_record = SeqRecord(seq)
        seq_record.id = '{0}:{1}-{2}'.format(seq_id, start+1, end)
        seq_record.description = ''
        seq_records.append(seq_record)
    return seq_records

def identify(genes, gene_db):
    prefix = generate_string()
    query_fn = '{}.query.fasta'.format(prefix)
    sbjct_fn = '{}.sbjct.fasta'.format(prefix)
    SeqIO.write(genes, query_fn, 'fasta')
    SeqIO.write(gene_db, sbjct_fn, 'fasta')
    result_fn = '{}.result.xml'.format(prefix)
    blastn_cline = NcbiblastnCommandline(
        query=query_fn, subject=sbjct_fn, evalue=1e-10,
        out=result_fn, outfmt=5, max_target_seqs=1,
    )
    stdout, stderr = blastn_cline()
    rows = []
    with open(result_fn, 'r') as result_file:
        blastn_records = NCBIXML.parse(result_file)
        for record in blastn_records:
            query_seqid, location = record.query.split(':')
            for aln in record.alignments:
                sbjct_seqid, description = aln.hit_def.split(' ', 1)
                sbjct_length = aln.length
                for hsp in aln.hsps:
                    ident_length = hsp.identities
                    align_length = hsp.align_length
                    if hsp.sbjct_start > hsp.sbjct_end:
                        query_seqid, location = record.query.split(':')
                        record.query = ':'.join([query_seqid, 'c'+location])
                    else:
                        pass
                    break
                break
            rows.append([sbjct_seqid, record.query, 
                        '{:.2f}'.format(ident_length / align_length * 100),
                        '{:.2f}'.format(align_length / sbjct_length * 100),
                        description])
    clean(prefix)
    return rows
        
def main():
    args = parse_args()
    genome = list(SeqIO.parse(args.input, 'fasta'))
    gene_db = list(SeqIO.parse(args.gene_db, 'fasta'))
    intervals = locate(genome, gene_db)
    output_writer = csv.writer(args.output, delimiter='\t')
    if intervals:
        genes = get_seq(genome, intervals)
        identify(genes, gene_db)        
        for row in identify(genes, gene_db):
            output_writer.writerow(row)

if __name__ == '__main__':
    sys.exit(main())
