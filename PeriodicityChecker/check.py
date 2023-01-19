'''
The main periodicity checker module. Calls the other modules to check for periodicity.

'''

import argparse
import os

from collapse_fastq_to_single_fasta import collapse
from build_contigs import build_contigs
from calculate_periodicity import calculate_periodicity


def main(args):
    '''
    Main function to bridge between argparse and the rest of the code.
    '''
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # collapse the fastq to a single fasta
    collapse_fastq_to_single_fasta(args.input_file, args.output_dir + '/collapsed.fa')

    # build contigs from the collapsed fasta
    contigs = build_contigs(args.output_dir + '/collapsed.fa', args.output_dir, 21)

    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Periodicity checker')
    parser.add_argument('input_file', type=str, help='Input fastq')
    parser.add_argument('output_dir', type=str, help='Output directory')
    args = parser.parse_args()