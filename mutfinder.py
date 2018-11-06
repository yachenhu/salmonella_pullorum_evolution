#!/usr/bin/python3

# Author: Yachen Hu (yachenhu@outlook.com)
# Dependencies:
#    Biopython	http://biopython.org/

"""
Detect site mutations and effects on product
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
from Bio.Blast.Applications import NcbiblastnCommandline, NcbiblastpCommandline
from Bio.Blast import NCBIXML
import csv

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gene_db', required=True, type=argparse.FileType('r'),
                        help="FASTA file of reference orfs")
    parser.add_argument('--input', required=True, type=argparse.FileType('r'),
                        help="FASTA file of genome sequence")
    parser.add_argument('--output', type=argparse.FileType('w'),
                        default=sys.stdout, help="Output file name")
    parser.add_argument('--threads', type=int, default=12, 
                        help="Number of threads to use in blast")
    return parser.parse_args()

def generate_string(size=8, chars=ascii_letters+digits):
    return ''.join(random.sample(chars, k=size))

def clean(prefix):
    for i in glob.glob('{0}.*'.format(prefix)):
        os.remove(i)

def find_allele(genome, orf):
    prefix = generate_string()
    query_fn = '{}.query.fasta'.format(prefix)
    sbjct_fn = '{}.sbjct.fasta'.format(prefix)
    SeqIO.write(genome, query_fn, 'fasta')
    SeqIO.write(orf, sbjct_fn, 'fasta')
    blastn_cline = NcbiblastnCommandline(
        query=query_fn, subject=sbjct_fn, evalue=1e-5,
        outfmt="'6 qseqid qstart qend sstart send'"
    )
    stdout, stderr = blastn_cline()
    csv_reader = csv.reader(stdout.splitlines(), delimiter='\t') # csv can handle any iterator
    qseqid, qstart, qend, sstart, send = next(csv_reader) # omit splits or duplicates
    qstart = int(qstart)
    qend = int(qend)
    sstart = int(sstart)
    send = int(send)
    location = '{0}-{1}'.format(qstart, qend)
    allele_seq = SeqIO.to_dict(genome)[qseqid][qstart-1:qend].seq
    if sstart > send:
        location = 'c' + location
        allele_seq = allele_seq.reverse_complement()
    else:
        pass
    allele = SeqRecord(allele_seq)
    allele.id = '{0}:{1}'.format(qseqid, location)
    allele.description = ''
    clean(prefix)
    return allele

def find_vars(query, sbjct):
    prefix = generate_string()
    query_fn = '{}.query.fasta'.format(prefix)
    sbjct_fn = '{}.sbjct.fasta'.format(prefix)
    SeqIO.write(query, query_fn, 'fasta')
    SeqIO.write(sbjct, sbjct_fn, 'fasta')
    result_fn = '{}.result.xml'.format(prefix)
    blastp_cline = NcbiblastpCommandline(
        query=query_fn, subject=sbjct_fn, evalue=1e-5, 
        out=result_fn, outfmt=5, 
    )
    stdout, stderr = blastp_cline()
    with open(result_fn, 'r') as result_file:
        record = NCBIXML.read(result_file)
        variants = []
        for aln in record.alignments:
            sbjct_length = aln.length
            for hsp in aln.hsps:
                ident_length = hsp.identities
                align_length = hsp.align_length
                for i, c in enumerate(hsp.match):
                    if c == ' ' or c == '+':
                        pos = i + 1
                        ref = hsp.sbjct[i]
                        alt = hsp.query[i]
                        variants.append('{0}{1}{2}'.format(ref, pos, alt))
                    else:
                        pass
                break # omit other hsps
            break # omit other alns
    clean(prefix)
    return ('{:.2f}'.format(ident_length / align_length * 100),
            '{:.2f}'.format(align_length / sbjct_length * 100), 
            ','.join(variants))

def main():
    args = parse_args()
    output_writer = csv.writer(args.output, delimiter='\t')
    genome = list(SeqIO.parse(args.input, 'fasta'))
    for orf in SeqIO.parse(args.gene_db, 'fasta'):
        allele = find_allele(genome, orf)
        if allele:
            query = SeqRecord(allele.seq.translate(to_stop=True))
            sbjct = SeqRecord(orf.seq.translate())
            identity, coverage, variants = find_vars(query, sbjct)
            output_writer.writerow([orf.id, allele.id, identity, 
                                    coverage, variants])
        else:
            pass

if __name__ == '__main__':
   sys.exit(main())
