import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
#from sourmash.fig import load_matrix_and_labels
import numpy as np
#import pandas as pd


def jaccard_to_evoldist(jaccard, ksize, b1=1.0, b2=1.0):
    # proportion of observed differences
    if jaccard ==0:
        return 1.0 # what should this be? (jaccard of 1.0 returns 0)
    p = 1 - np.power(2*jaccard/(jaccard + 1),(1/float(ksize)))
    # corrected evolutionary distance
    d = -(b1*np.log((1-p)/b2))
    return d


def main(args):
    # from query csv or np matrix, build dictionary of group:: filenames
    jaccard_dists = np.load(open(args.jaccard_np, 'rb'))
    jaccard_to_evoldist_vectorized = np.vectorize(jaccard_to_evoldist)
    # convert each jaccard dist value, ignoring 1.0, 0.0 values
    evoldists = jaccard_to_evoldist_vectorized(jaccard_dists, args.ksize)
    with open(args.evoldist_np, 'wb') as fp:
        np.save(fp, evoldists)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("jaccard_np")
    p.add_argument("--ksize", type=int, required=True)
    p.add_argument("--evoldist_np", required=True)
    p.add_argument("--b1", type=float, default=1.0)
    p.add_argument("--b2", type=float, default=1.0)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
