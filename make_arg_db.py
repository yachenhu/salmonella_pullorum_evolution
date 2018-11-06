#!/usr/bin/python3

# Author: Yachen Hu (yachenhu@outlook.com)
# Dependencies:
#    Biopython  http://biopython.org/


"""
Create ARG database by using resfinder files
"""

import argparse
import sys
import os
import csv
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--note', required=True, type=argparse.FileType('r'),
                        help="Resfinder notes file name")
    parser.add_argument('--output', type=argparse.FileType('w'), 
                        default=sys.stdout, help="Output sequence file")
    parser.add_argument('input', nargs='+', help="Input sequence files")
    return parser.parse_args()

def get_res_dict(res_file):
    csv_reader = csv.reader(res_file, delimiter=':')
    res_dict = {}
    for row in csv_reader:
        if not row[0].startswith('#'):
            res_dict[row[0]] = row[1]
        else:
            pass
    return res_dict

def main():
    args = parse_args()
    res_dict = get_res_dict(args.note)
    for gene_fn in args.input:
        for gene in SeqIO.parse(gene_fn, 'fasta'):
            name = gene.id.split('_', 1)[0]
            gene.description = res_dict[name]
            SeqIO.write(gene, args.output, 'fasta')

if __name__ == '__main__':
    sys.exit(main())
