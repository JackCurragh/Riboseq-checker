'''
read in tsv file and calculate the periodicity of each contig

input is the path to a tsv file with the following columns:
    - contig
    - start
    - end
    - count
'''

import pandas as pd
import argparse



def calculate(bed_path:str, output_path:str):
    # read in tsv and use the following list as columns rather than the first row
    df = pd.read_csv(bed_path, sep='\t', header=None, names=['contig', 'start', 'end', 'count'])
    contigs = df['contig'].unique()
    with open(output_path, 'w') as f:

        # Initialize the output file
        f.write('# Periodicity Checker')
        f.write('# oeriodicity is calculated as the number of reads in the most abundant frame divided by the total number of reads')
        f.write('contig\tframe0\tframe1\tframe2\ttotal\tperiodicity')

        # Iterate through each contig and calculate the periodicity
        for contig in contigs:
            frame_counts = {0:0, 1:0, 2:0}
            contig_df = df[df['contig'] == contig]
            for i in range(1, max(contig_df['start'])):
                if i in contig_df['start'].values:
                    frame_counts[i%3] += contig_df[contig_df['start'] == i]['count'].values[0]
            score = max(frame_counts.values()) / sum(frame_counts.values())
            f.write(f'{contig}\t{frame_counts[0]}\t{frame_counts[1]}\t{frame_counts[2]}\t{frame_counts[0]+frame_counts[1]+frame_counts[2]}\t{score}')
    f.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate the periodicity of each contig in a tsv file.')
    parser.add_argument('--tsv', type=str, help='Path to the tsv file.')
    parser.add_argument('--output', type=str, help='Path to the output file.')
    args = parser.parse_args()

    calculate(args)