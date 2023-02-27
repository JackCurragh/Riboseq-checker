"""
Obtain a score for the uniformity of a ribosome profile
"""

import argparse
import pandas as pd
from scipy.stats import kstest, uniform
import numpy as np
import math

def main(args):
    """
    wrapper function to obtain the uniformity score
    """
    # Read in the bed file
    df = pd.read_csv(args.bed, sep='\t', header=None)
    print(df)
    # Get the read counts
    for transcript in df[0].unique():
        print(transcript)
        # Get the read counts for this transcript
        df_transcript = df[df[0] == transcript]
        # Get the read counts
        counts = np.array(df_transcript[3].values)
        # calculate total number of reads
        total_reads = sum(counts)
        # calculate probability distribution of coverage values
        probability_distribution = [float(count) / total_reads for count in counts]

        max_entropy = -sum(p * math.log(p, 2) for p in probability_distribution if p != 0)
        print(max_entropy)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Obtain a score for the uniformity of a ribosome profile')
    parser.add_argument('--bed', help='bed file containing read counts')
    args = parser.parse_args()
    main(args)
