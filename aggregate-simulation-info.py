import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
import numpy as np
import pandas as pd


def process_infofile(inF):
    # Read and modify single simulation csv
    seq_basename = os.path.basename(inF).rsplit(".tsv.xz")[0]
    info_csv = pd.read_csv(inF, sep = "\t", compression="xz")
    info_csv["name"] = seq_basename + "-seed" + info_csv["seed"].apply(str)
    info_csv.drop(columns=["seq1", "seq2"], inplace=True)
    info_csv.set_index("name", inplace=True)
    return info_csv


def main(args):
    # get basename for these sequences
    filelist = [x.strip() for x in open(args.info_filelist, "r")]
    frames = [process_infofile(f) for f in filelist]
    full_info = pd.concat(frames)
    full_info.to_csv(args.output_csv, index_label="name")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("info_filelist")
    p.add_argument("output_csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
