import os
import sys
import argparse
import glob
import pprint

import screed
import sourmash
import numpy as np
import pandas as pd

NUCLEOTIDE_KSIZES=[21,31,51]
PROTEIN_KSIZES=[8,9,10,11,12]
DAYHOFF_KSIZES=[15,16,17,18,19]
HP_KSIZES=[33,35,37,39,41,42]
NUCLEOTIDE_SCALED = 1000
PROTEIN_SCALED = 100
DAYHOFF_SCALED = 100
HP_SCALED = 100

def determine_appropriate_fresh_minhash(alphabet, ksize, scaled_val, ignore_abundance=False):
    # default behavior is to track abundance
    abund = not ignore_abundance
    if alphabet == "nucleotide":
        mh = sourmash.MinHash(ksize=ksize, n=0, scaled=scaled_val, track_abundance=abund, is_protein=False)
    elif alphabet == "protein":
        k=ksize*3 ## need to multiply bt 3 to get same ksize, bc add_protein method does k/3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=False)
    elif alphabet == "dayhoff":
        k=ksize*3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=True, hp=False)
    elif alphabet == "hp":
        k=ksize*3
        mh = sourmash.MinHash(ksize=k, n=0, scaled=scaled_val, track_abundance=abund, is_protein=True, dayhoff=False, hp=True)
    return mh


def build_sig_from_file(input_file, alphabet, ksize, scaled, ignore_abundance=False, translate=False):
    sig=""
    records = records=screed.open(input_file)
    signame = os.path.basename(input_file).rsplit(".fasta")[0]
    # start with fresh minhash
    mh = determine_appropriate_fresh_minhash(alphabet, ksize, scaled, ignore_abundance)
    if records:
        for record in records:
            if alphabet == "nucleotide" or translate:
                mh.add_sequence(record.sequence, force=True)
            else:
                mh.add_protein(record.sequence)
        # minhash --> signature
        sig = sourmash.SourmashSignature(mh, name=signame)
    return sig


def convert_to_pdist(corr, ksize):
    # return proportion of observed differences
    if corr ==0:
        return 1.0 # what should this be? (jaccard of 1.0 returns 0)
    p = 1 - np.power(2*corr/(corr + 1),(1/float(ksize)))
    return p


def compare_sequences(row, fasta_dir, input_alpha, alphabets, ksizes, scaled_vals):
    seq1_fastaF = os.path.join(fasta_dir, row["name"] + "-seq1" + ".fasta")
    seq2_fastaF = os.path.join(fasta_dir, row["name"] + "-seq2" + ".fasta")
    translate = False
    for alpha in alphabets:
        if input_alpha == "nucleotide" and alpha != "nucleotide":
            translate=True
        for ksize in ksizes[alpha]:
            for scaled in scaled_vals[alpha]:
                jaccard_name = f"{alpha}-k{str(ksize)}-scaled{str(scaled)}.jaccard"
                jaccard_pdist = f"{alpha}-k{str(ksize)}-scaled{str(scaled)}.jaccard-pdist"
                containment_name = f"{alpha}-k{str(ksize)}-scaled{str(scaled)}.containment"
                containment_pdist = f"{alpha}-k{str(ksize)}-scaled{str(scaled)}.containment-pdist"
                # build signatures
                seq1_sig = build_sig_from_file(seq1_fastaF, alpha, ksize, scaled, ignore_abundance=True, translate=translate)
                seq2_sig = build_sig_from_file(seq2_fastaF, alpha, ksize, scaled, ignore_abundance=True, translate=translate)
                # compare (jaccard, containment)
                jaccard = seq1_sig.jaccard(seq2_sig)
                contain1 = seq1_sig.contained_by(seq2_sig)
                contain2 = seq2_sig.contained_by(seq1_sig)
                max_contain = max(contain1,contain2)
                # convert to pdists
                j_pdist = convert_to_pdist(jaccard, ksize)
                c_pdist = convert_to_pdist(max_contain, ksize)
                row[jaccard_name] = jaccard
                row[jaccard_pdist] = j_pdist
                row[containment_name] = max_contain
                row[containment_pdist] = c_pdist
    return row


def process_infofile(inF):
    # Read and modify single simulation csv
    seq_basename = os.path.basename(inF).rsplit(".tsv")[0]
    info_csv = pd.read_csv(inF, sep = "\t", compression="xz")
    info_csv["name"] = seq_basename + "-seed" + info_csv["seed"].apply(str)
    info_csv.drop(columns=["seq1", "seq2"], inplace=True)
    return info_csv

def main(args):
    # get basename for these sequences
    info_csv = process_infofile(args.simulation_csv)
    alphabets_to_compare = ["nucleotide", "protein", "dayhoff", "hp"]
    ksize_info = {"nucleotide": [21,31,51], "protein": [7,8,9,10,11,12], "dayhoff": [15,16,17,18,19], "hp": [33,35,37,39,42]}
    scaled_info = {"nucleotide": [1000], "protein": [100], "dayhoff": [100], "hp": [100]}
    info_csv.apply(compare_sequences, axis=1, args=(str(args.fasta_dir), str(args.fasta_alphabet), alphabets_to_compare, ksize_info, scaled_info))
    #write to csv, don't write index
    info_csv.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--simulation-csv")
    p.add_argument("--fasta-dir")
    p.add_argument("--fasta-alphabet", default="nucleotide")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

