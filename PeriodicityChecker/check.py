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
   
    print("### Periodicity Checker ###")
    print(f"Checking the periodicity of the following file: {args.i}")
    print(f"Output will be written to: {args.output}")

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # collapse the fastq to a single fasta
    fa_path = collapse(args.i, args.output + '/collapsed.fa')

    # build contigs from the collapsed fasta
    contigs = build_contigs(fa_path, args.output)

    # align reads back to contigs
    bam = realign(contigs, fa_path, args.output)

    # convert bam to bed
    bed = generate_profile(bam, contigs, offset=15)

    # calculate periodicity
    periodicity = calculate(bed, contigs)



    pass

if __name__ == "__main__":
    '''
    This tool has 3 modes:
    1. Check the periodicity of a single fastq file without knowing the organism of origin
    2. Check the periodicity of a single fastq file knowing the organism of origin
    3. Check the periodicity of a BAM file (already aligned)
    '''

    # parser = argparse.ArgumentParser(description='Periodicity checker')
    # parser.add_argument('-A', type=str, help='Mode A: "Reference Agnostic" - Check the periodicity of a single fastq file without knowing the organism of origin')
    # parser.add_argument('-O', type=str, help='Mode O: "Organism Known" - Check the periodicity of a single fastq file knowing the organism of origin')

    # parser.add_argument('-i', type=str, help='Input fastq')
    # parser.add_argument('-o', type=str, help='Output directory')

    # parser.add_argument('-v', type=str, default=False, help='Verbose')
    # args = parser.parse_args()
    # main(args)


    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-A", "--agnostic", action="store_true", help='Mode A: "Agnostic" - Check the periodicity of a single fastq file without knowing the organism of origin')
    group.add_argument("-O", "--organism", action="store_true", help='Mode O: "Organism Known" - Check the periodicity of a single fastq file knowing the organism of origin')
    group.add_argument("-B", "--BAM", action="store_true", help='Mode B: "BAM input" Check the periodicity of a BAM file (already aligned)')
    parser.add_argument("-q", "--fastq", help="fastq name")
    parser.add_argument("-g", "--gff", help="gff file name")
    parser.add_argument("-b", "--bam", help="bam file name")
    parser.add_argument('--output', type=str, help='Output directory')

    args = parser.parse_args()

    if args.A and not args.q:
        parser.error("-A mode requires -f option.")
    if args.O and (not args.q or not args.g):
        parser.error("-O mode requires -f and -g options.")
    if args.B and not args.b:
        parser.error("-B mode requires -b option.")

    main(args)