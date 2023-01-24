'''
Script that takes a bed file and calculates the fourier transform of the read counts


'''

import pandas as pd
import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt


def fourier(bed_path):
    '''
    Calculate the fourier transform of the read counts in a bed file
    '''
    # Read in the bed file
    df = pd.read_csv(bed_path, sep="\t", header=None, names=['contig', 'start', 'end', 'count'])

    # Extract the start positions
    grouped = df.groupby('contig')
    frame_counts = {0:0, 1:0, 2:0}

    for name, group in grouped:
        for i in range(1, max(group['start'])):
            if i in group['start'].values:
                frame_counts[i%3] += group[group['start'] == i]['count'].values[0]

        plt.plot(group['start'], group['count'], label=name)
        print(frame_counts)
        plt.xlabel('Position')
        plt.ylabel('Counts')
        plt.legend()
        plt.show()
    # for contig in df['contig'].unique():
    #     contig_df = df[df['contig'] == contig]
    #     counts = contig_df['count'].tolist()

    #     plt.plot(counts)
    #     plt.xlabel("position")
    #     plt.ylabel("coutns")
    #     plt.show()
        # # Perform the FFT
        # fft_result = fft(counts)

        # # Get the frequency amplitudes
        # frequency_amplitudes = np.abs(fft_result)

        # # Plot the frequency amplitudes
        # plt.plot(frequency_amplitudes)
        # plt.xlabel("Frequency")
        # plt.ylabel("Amplitude")
        # plt.show()

if __name__ == "__main__":
    fourier('example/SRR2064438.sorted.bam.bed.sorted')