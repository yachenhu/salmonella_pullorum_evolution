#!/usr/bin/python3

"""
Detect deletion structure variations (DSV) by comparing two sequences
"""

import argparse
import sys
import os
import tempfile
from subprocess import run, PIPE, CalledProcessError
import csv

from Bio import SeqIO
from pybedtools import BedTool, Interval, create_interval_from_list

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '-query', required=True, type=str, metavar='FILE',
        help="Query FASTA filename"
    )
    parser.add_argument(
        '-reference', required=True, type=str, metavar='FILE',
        help="Reference FASTA filename"
        )
    parser.add_argument(
        '-length', default=75, type=float, metavar='FLOAT',
        help="Minimum length of match [75]"
        )
    parser.add_argument(
        '-identity', default=90, type=float, metavar='FLOAT',
        help="Minumu identity of match [90]"
        )
    parser.add_argument(
        '-out', default=sys.stdout, type=argparse.FileType('w'), metavar='FILE',
        help="output gap BED filename"
        )
    return parser.parse_args()

def nucmer(ref_filename, qry_filename):
    with tempfile.TemporaryDirectory() as tmpdirname:
        nucmer_cmd = ['nucmer', '--mum', '-p', os.path.join(tmpdirname,'out'),
                      ref_filename, qry_filename]
        nucmer_proc = run(nucmer_cmd, stdout=PIPE, stderr=PIPE, check=True)
        show_coords_cmd = ['show-coords', '-H', '-T', '-r',
                           os.path.join(tmpdirname, 'out.delta')]
        show_coords_proc = run(show_coords_cmd, check=True,
                               stdout=PIPE, stderr=PIPE)
        return(show_coords_proc.stdout.decode())

def filter_coords(coords, length, identity):
    """Filter NUCmer coords."""
    for coord in coords:
        if float(coord[4]) >= length and float(coord[6]) >= identity:
            yield coord

def coords_to_bed(coords):
    """Convert NUCmer coords to BedTool object and merge."""
    intervals = []
    for row in coords:
        interval = create_interval_from_list(
            [row[7], int(row[0])-1, int(row[1])]
            )
        intervals.append(interval)
    return BedTool(intervals).sort().merge()

def genome_to_bed(genome):
    """Generate genome BedTool for genome sequence file."""
    intervals = []
    for seqrecord in SeqIO.parse(genome, 'fasta'):
        interval = create_interval_from_list(
            [seqrecord.id, 0, len(seqrecord)]
            )
        intervals.append(interval)
    return BedTool(intervals)

def main():
    args = parse_args()
    try:
        nucmer_out = nucmer(args.reference, args.query)
        nucmer_coords = csv.reader(nucmer_out.splitlines(), delimiter='\t')
        valid_coords = filter_coords(nucmer_coords, args.length, args.identity)
        match_bed = coords_to_bed(valid_coords)
        genome_bed = genome_to_bed(args.reference)
        deletion_bed = genome_bed.subtract(match_bed)
    except CalledProcessError as e:
        return e
    writer = csv.writer(args.out, delimiter='\t', lineterminator='\n')
    for row in deletion_bed:
        writer.writerow(row)

if __name__ == '__main__':
    sys.exit(main())
