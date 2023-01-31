import argparse
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator


def collapse(fastq: str, fasta: str) -> str:
    '''
    Collapse a FASTQ file to a FASTA file with read counts in the header.
    '''
    # Store the unique reads in a dictionary with their counts
    unique_reads = {}

    # Open the FASTQ file for reading, handling gzip if necessary
    # by checking the magic number at the beginning of the file
    with open(fastq, 'rb') as f:
        magic_number = f.read(2)
        if magic_number == b'\x1f\x8b':
            f = gzip.open(fastq, 'rt')
        else:
            f.seek(0)
            f = open(fastq, 'r')

        for title, sequence, quality in FastqGeneralIterator(f):
            if sequence in unique_reads:
                unique_reads[sequence]['count'] += 1
            else:
                unique_reads[sequence] = {}
                unique_reads[sequence]['count'] = 1

    # Open the FASTA file for writing, handling gzip if necessary
    if fasta.endswith('.gz'):
        f = gzip.open(fasta, 'wt')
    else:
        f = open(fasta, 'w')

    # Write the unique reads to the FASTA file
    read_number = 1
    for seq, vals in unique_reads.items():
        f.write(f'>read{read_number}_x{vals["count"]}\n')
        f.write(f"{seq}\n")
        read_number += 1
    return fasta


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', help='path to input FASTQ file')
    parser.add_argument('-o', help='path to output collapsed FASTA file (add .gz to the end to compress)')

    args = parser.parse_args()
    collapse(args.i, args.o)