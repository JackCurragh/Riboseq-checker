'''
This script specifically handles the bam file input. It is called by the main script: check.py

Frame bias is calculated by the alignments in the bam file alone
'''


import pysam
import subprocess
import os



def is_name_sorted(bam_path) -> bool:
    '''
    check if a bam file is name sorted
    '''
    bam = pysam.Samfile(bam_path, 'rb')
    prev_read = None
    for read in bam:
        if prev_read:
            if read.qname < prev_read.qname:
                return False
        prev_read = read
    return True


def sort_by_name(bam_path, num_threads=1) -> str:
    '''
    sort a bam file by name using samtools
    '''
    sorted_bam_path = bam_path.replace(".bam", "_sorted.bam")
    subprocess.call(["samtools", "sort", "-@", str(num_threads), "-n", bam_path, '-o', sorted_bam_path])
    return sorted_bam_path


def is_indexed(bam_path) -> bool:
    '''
    check if a bam.bai file exists
    '''
    return os.path.exists(bam_path + ".bai")


def index_bam(bam_path) -> None:
    '''
    index a bam file using samtools
    '''
    subprocess.call(["samtools", "index", bam_path])


def check_bam(bam_path, output_path, num_threads=1, offset=15) -> str:
    '''
    calculate frame bias from a bam file

    inputs:
        bam_path: path to bam file
        num_threads: number of threads to use

    outputs:
        frame_counts: dictionary with frame bias for each transcript/reference
    '''
    # check if bam file is sorted by name
    # if not is_name_sorted(bam_path):
    #     print("bam file not sorted by name, sorting now...")
    #     bam_path = sort_by_name(bam_path, num_threads)

    if not is_indexed(bam_path):
        print("bam file not indexed, indexing now...")
        bam_path = index_bam(bam_path)

    with open(output_path, 'w') as f:

        # Initialize the output file
        f.write('# Periodicity Checker\n')
        f.write('# periodicity is calculated as the number of reads in the most abundant frame divided by the total number of reads\n')
        f.write('contig\tframe0\tframe1\tframe2\ttotal\tperiodicity\n')

        bam = pysam.Samfile(bam_path, 'rb')

        frame_counts = {}

        # get frame bias per transcript
        for transcript in bam.references:
            frame_counts[transcript] = {0:0, 1:0, 2:0}
            transcript_reads = bam.fetch(transcript)
            for read in transcript_reads:
                if read.is_unmapped:
                    continue
                frame_counts[transcript][(read.reference_start + offset)%3] += int(read.qname.split('_x')[1])

            score = max(frame_counts[transcript].values()) / sum(frame_counts[transcript].values())
            f.write(f'{transcript}\t{frame_counts[transcript][0]}\t{frame_counts[transcript][1]}\t{frame_counts[transcript][2]}\t{frame_counts[transcript][0]+frame_counts[transcript][1]+frame_counts[transcript][2]}\t{score}\n')
    f.close()   

    return frame_counts