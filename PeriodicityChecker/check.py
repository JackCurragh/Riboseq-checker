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
from bam_check import check_bam


def run_agnostic(args, fa_path):
    '''
    wrapper function for agnostic mode. This invloves constructing contigs from the 
    reads with inchworm and then aligning the reads back to the contigs with bowtie.
    periodicity can be obtained from the resulting profilies

    Necessary arguments:
        -q: path to fastq file (but this is handled prior to this function)
        --output: path to output directory
    '''
    # build contigs from the collapsed fasta
    contigs = build_contigs(fa_path, args.output)

    # align reads back to contigs
    bam = realign(contigs, fa_path, args.output)

    # convert bam to bed
    bed = generate_profile(bam, contigs, offset=15)

    # calculate periodicity
    periodicity = calculate(bed, contigs)
    return periodicity


def run_organism(args, fa_path):
    '''
    wrapper function for organism mode. This involves aligning the reads to a 
    reference with bowtie and then calculating periodicity from the resulting
    profiles

    Necessary arguments:
        -f or --fasta: path to reference fasta file
        --output: path to output directory
        -q: path to fastq file (but this is handled prior to this function)

    '''
    # get the path to the reference
    if args.f:
        reference_path = args.f
    else:
        reference_path = args.fasta

    # align reads to reference
    bam = realign(reference_path, fa_path, args.output)

    # convert bam to bed
    bed = generate_profile(bam, reference_path, offset=15) 

    # calculate periodicity
    periodicity = calculate(bed, reference_path)
    return periodicity


def run_bam(args):
    '''
    Wrapper function for bam mode. This involves calculating periodicity from a
    pre-existing bam file
    '''
    check_bam(args.bam, f"{args.output}/periodicity.txt")



def main(args):
    '''
    Main function to bridge between argparse and the rest of the code.
    '''
   
    print("### Periodicity Checker ###")
    print(f"Output will be written to: {args.output}")

    if not os.path.exists(args.output):
        os.makedirs(args.output)


    if args.agnostic:
        collapsed_fa_path = collapse(args.i, args.output + '/collapsed.fa')
        periodicity = run_agnostic(args, collapsed_fa_path)

    elif args.organism:
        collapsed_fa_path = collapse(args.i, args.output + '/collapsed.fa')
        periodicity = run_organism(args, collapsed_fa_path)

    elif args.bam:
        periodicity = run_bam(args)

    else:
        raise Exception("No mode selected")




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
    parser.add_argument("-q", "--fastq", help="fastq path")
    parser.add_argument("-f", "--fasta", help="reference fasta file path")
    parser.add_argument("-b", "--bam", help="bam file path")
    parser.add_argument('--output', type=str, help='Output directory')

    args = parser.parse_args()
    print(args)

    if args.agnostic and not args.q:
        parser.error("-A mode requires -f option.")
    if args.organism and (not args.q or not args.g or not args.f):
        parser.error("-O mode requires -f, -q and -g options.")
    if args.BAM and not args.bam:
        parser.error("-B mode requires -b option.")

    main(args)