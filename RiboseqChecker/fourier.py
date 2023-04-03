'''
Script that takes a bed file and calculates the fourier transform of the read counts


'''
import argparse
import pandas as pd
import numpy as np
from scipy.fftpack import fft, fftfreq
import matplotlib.pyplot as plt

def plot_profile(name, frame_counts):
    '''
    plot the ribosome profile given read counts
    '''
    for i, frame in enumerate(frame_counts):
        plt.plot(frame_counts[frame], label=f'Frame {i}')

    score = max(sum(frame_counts[frame]) for frame in frame_counts) / sum(sum(frame_counts[frame]) for frame in frame_counts)
    tx_name = name.split('|')[0]
    plt.title(f'{tx_name}\n{score}')
    plt.xlabel('Position')
    plt.ylabel('Counts')
    plt.legend()
    plt.show()


def plot_fourier(fftfreq_result, frequency_amplitudes):
    '''
    Visualise a fourier transform of a signal using matplotlib
    '''
    plt.plot(fftfreq_result, frequency_amplitudes, marker='o')
    plt.xlabel("Frequency")
    plt.ylabel("Amplitude")
    plt.show()


def fourier(counts):
    '''
    Calculate the fourier transform of the read counts in a bed file
    '''

    # # Perform the FFT
    fft_result = fft(counts)
    fftfreq_result = fftfreq(len(counts), 1)

    # Get the frequency amplitudes
    frequency_amplitudes = np.abs(fft_result)

    pos_indexes = [i for i, val in enumerate(fftfreq_result) if val > 0.05]
    pos_result = frequency_amplitudes[pos_indexes]

    for x, y in zip(fftfreq_result, frequency_amplitudes):
        if y == max(pos_result) and y == min(pos_result):

            print('Warning: Periodicity score is 0', y, max(pos_result), min(pos_result))
        if y == max(pos_result):
            if x > 0:
                if y != min(pos_result):
                    index = list(fftfreq_result).index(x)
                    full_list = frequency_amplitudes[index-10: index+10]
                    surrounding_list_below = frequency_amplitudes[index-10: index-3]
                    surrounding_list_above = frequency_amplitudes[index+3: index+10]
                    try:
                        surrounding_list = surrounding_list_below + surrounding_list_above
                        if x > 0.32 and x < 0.34:
                            return (max(full_list) - (sum(surrounding_list)/len(surrounding_list)))/max(full_list)  
                        else:
                            return 0
                    except:
                        plot_fourier(fftfreq_result, frequency_amplitudes)
                        raise ValueError('Error in calculating periodicity score')
                else:
                    raise ValueError('Y also equals min')

    return 0

def sort_dict_by_value(d):
    '''
    Sort a dictionary by value
    '''
    return {k: v for k, v in sorted(d.items(), key=lambda item: item[1])}

def main(args):
    '''
    Main function to bridge between argparse and the rest of the code.
    '''
    df = pd.read_csv(args.bed, sep="\t", header=None, names=['contig', 'start', 'end', 'count'])

    # Extract the start positions
    grouped = df.groupby('contig')
    #calculate totals per contig from grouped dataframe
    totals_per_contig = {}
    for name, group in grouped:
        if 'protein_coding' in name:
            totals_per_contig[name] = sum(group['count'])
    
    top_contigs = list(sort_dict_by_value(totals_per_contig).keys())[-100:]
    periodicity_scores = {}
    for name, group in grouped:
        if name in top_contigs:
            counts = []
            frame_counts = {0:[], 1:[], 2:[]}

            for i in range(1, max(group['start'])):
                if i in group['start'].values:
                    frame_counts[i%3].append(group[group['start'] == i]['count'].values[0])
                    counts.append(group[group['start'] == i]['count'].values[0])
                else:
                    frame_counts[i%3].append(0)
                    counts.append(0)

            score = fourier(counts)
            if score > 0:
                periodicity_scores[name] = score
        
    print(sum( [periodicity_scores[x] for x in periodicity_scores] ) / len(periodicity_scores))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate periodicity from a bed file')
    parser.add_argument('--bed', help='Path to bed file')
    args = parser.parse_args()
    main(args)


'''
 python PeriodicityChecker/fourier.py --bed  example/SRR2064438.sorted.bam.bed.sorted
 
 '''