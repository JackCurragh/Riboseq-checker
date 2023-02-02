"""
Obtain a score for the uniformity of a ribosome profile
"""

import argparse
import pandas as pd
from scipy.stats import kstest, uniform
import numpy as np

def main(args):
    """
    wrapper function to obtain the uniformity score
    """
    rng = np.random.default_rng()
    print(kstest(uniform.rvs(size=5, random_state=rng), [1,2,3,4,5]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Obtain a score for the uniformity of a ribosome profile')
    parser.add_argument('--bed', help='bed file containing read counts')
    args = parser.parse_args()
    main(args)
