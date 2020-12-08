import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
import numpy as np
import pandas as pd


def seq_to_file(row, outdir):
    seq1 = row["seq1"].replace("-", "")
    seq2 = row["seq2"].replace("-", "")
    basename = row["name"]
    name1 = basename + "_seq1"
    file1 = os.path.join(outdir, name1 + ".fasta")
    name2 = basename + "_seq2"
    file2 = os.path.join(outdir, name2 + ".fasta")
    with open(file1, "w") as f1:
        f1.write(">" + name1 + "\n" + seq1 + "\n")
    with open(file2, "w") as f2:
        f2.write(">" + name2 + "\n" + seq2 + "\n")


def main(args):
    # get basename for these sequences
    seq_basename = os.path.basename((args.input_tsv)).rsplit(".tsv")[0]
    # Read simulation csv
    info_csv = pd.read_csv(args.input_tsv, sep = "\t", compression="xz")
    #make name for each sequence
    info_csv["name"] = seq_basename + "-seed" + info_csv["seed"].apply(str)
    info_csv.apply(seq_to_file, axis=1, args=(str(args.outdir),))


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("input_tsv")
    p.add_argument("--outdir", default = "")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
