'''
The main periodicity checker module. Calls the other modules to check for periodicity.

'''

import argparse
import os

from collapse_fastq_to_single_fasta import collapse
from build_contigs import build_contigs
from realign_to_contigs import realign
from calculate_periodicity import calculate
from bam_to_ribosome_profile import generate_profile


def main(args):
    '''
    Main function to bridge between argparse and the rest of the code.
    '''
    if not os.path.exists(args.o):
        os.makedirs(args.o)

    # collapse the fastq to a single fasta
    fa_path = collapse(args.i, args.o + '/collapsed.fa')

    # build contigs from the collapsed fasta
    contigs = build_contigs(fa_path, args.o)

    # align reads back to contigs
    bam = realign(contigs, fa_path, args.o)

    # convert bam to bed
    bed = generate_profile(bam, contigs, offset=15)

    # calculate periodicity
    periodicity = calculate(bed, contigs)



    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Periodicity checker')
    parser.add_argument('-i', type=str, help='Input fastq')
    parser.add_argument('-o', type=str, help='Output directory')
    args = parser.parse_args()