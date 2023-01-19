'''
construct bowtie index from fasta file and realign reads to contigs
convert and sort alignments to bam file

inputs:
    contigs: path to contigs fasta file
    reads: path to reads fastq file
    tmp_dir: path to tmp directory
    
outputs:
    bam: path to bam file

'''
import subprocess
import argparse


def index_contigs(contigs: str, tmp_dir: str) -> str:
    '''
    construct bowtie index from fasta file and realign reads to contigs
    '''
    # construct bowtie index
    cmd = ["bowtie-build",  contigs, tmp_dir + '/contig_index']
    subprocess.call(cmd)


def algin_reads_to_contigs(reads: str, tmp_dir: str, num_threads=1) -> str:
    '''
    realign reads to contigs
    '''
    # realign reads to contigs
    cmd = [f"bowtie", '-a', '--norc', '-p', num_threads,  '-v', 3, '--seedlen', 25, tmp_dir + '/contig_index', '-q', reads, '-S', tmp_dir + '/contigs.sam']
    subprocess.call(cmd)

    cmd = [f"samtools", 'sort', '-@', num_threads, tmp_dir + '/contigs.sam', '-o', tmp_dir + '/contigs.bam']
    subprocess.call(cmd)

    return tmp_dir + '/contigs.bam'


def realign(contigs: str, reads: str, tmp_dir: str, num_threads=1) -> str:
    '''
    wrapper function for realigning reads to contigs

    inputs:
        contigs: path to contigs fasta file
        reads: path to reads fastq file
        tmp_dir: path to tmp directory
        num_threads: number of threads to use

    outputs:    
        bam: path to sorted bam file
    '''
    index_contigs(contigs, reads, tmp_dir)
    bam = algin_reads_to_contigs(reads, tmp_dir, num_threads)
    return bam
