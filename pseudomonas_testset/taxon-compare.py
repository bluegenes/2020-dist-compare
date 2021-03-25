import os
import sys
import argparse
import glob
import pprint

import pandas as pd

import screed
import sourmash
from sourmash.sourmash_args import load_file_as_signatures

from collections import defaultdict, namedtuple

CompareResult = namedtuple('CompareResult',
                           'comparison_name, anchor_name, ref_name, lowest_common_rank, anchor_sciname, alphabet, ksize, scaled, jaccard, max_containment, anchor_containment, anchor_hashes, query_hashes, num_common')

def compare_sigs(sigA, sigB, comparison_name, lowest_common_rank, anchor_sciname, alpha, ksize, scaled):
    sigA_numhashes = len(sigA.minhash.hashes)
    sigB_numhashes = len(sigB.minhash.hashes)
    intersect_numhashes = sigA.minhash.count_common(sigB.minhash)
    jaccard = sigA.jaccard(sigB)
    containA = sigA.contained_by(sigB)
    max_contain = sigA.max_containment(sigB)
    #max_contain = max(containA,containB)
    return CompareResult(comparison_name, str(sigA).split(" ")[0], str(sigB).split(" ")[0], lowest_common_rank, anchor_sciname, alpha, ksize, scaled, jaccard, max_contain, containA, sigA_numhashes, sigB_numhashes, intersect_numhashes)

def main(args):
    ksize=args.ksize
    scaled=args.scaled
    alphabet=args.alphabet
    sigext = args.sig_extension
    if alphabet == "nucleotide":
        moltype = "DNA"
    else:
        moltype = alphabet

    # find all sigs
    siglist = [x.rstrip() for x in open(args.siglist)]
    sigD={}
    for sigF in siglist:
        name = os.path.basename(sigF).rsplit(sigext)[0]
        if not os.path.exists(sigF):
            full_sigF = os.path.join(args.sigdir, sigF)
            if not os.path.exists(full_sigF):
                print(f"sig {name} cannot be found at {sigF} or within sigdir {args.sigdir}")
                continue
            else:
                sigF=full_sigF
        sigD[name] = sigF


    taxon_comparisons = []
    compareInfo = pd.read_csv(args.comparison_csv)
    # from comparison csv, build dictionary of group:: anchor, comparison filenames
    compareInfo = compareInfo.groupby(["anchor_acc", "lowest_common_rank", "anchor_sciname"]).agg({"compare_accs":lambda x: list(x)})#.reset_index()
    # loop through comparisons
    for n, (anchor_acc, lowest_common_rank, anchor_sciname) in enumerate(compareInfo.index):
        if n !=0 and n % 50 == 0:
            print(f"... assessing {n}th taxon comparison, lowest common rank: {lowest_common_rank}, anchor: {anchor_acc} {anchor_sciname}\n")

        # select and load anchor sig
        selector = load_file_as_signatures(sigD[anchor_acc], ksize=ksize, select_moltype=moltype)
        anchor_sig = next(selector)

        # iterate through comparison sigs
        compare_accs = compareInfo.loc[(anchor_acc, lowest_common_rank)]["compare_accs"].values[0]
        for compare_acc in compare_accs:
            # select and load comparison sig
            selector = load_file_as_signatures(sigD[compare_acc], ksize=ksize, select_moltype=moltype)
            compare_sig = next(selector)
            comparison = compare_sigs(anchor_sig, compare_sig, f"{anchor_acc}_x_{compare_acc}", lowest_common_rank, anchor_sciname, alphabet, ksize, scaled)
            taxon_comparisons.append(comparison)

    # convert path comparison info to pandas dataframe
    comparisonDF = pd.DataFrame.from_records(taxon_comparisons, columns = CompareResult._fields)

    # print to csv
    comparisonDF.to_csv(args.output_csv, index=False)
    print(f"done! taxon comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--comparison-csv", default="gtdb-r95-reps.pathinfo.tsv")
    p.add_argument("--lineages-csv", default="gtdb-r95-reps.lineages.protein-filenames.reordered.csv")
    p.add_argument("--siglist", default="gtdb95-evolpaths/gtdb95-evolpaths.signatures.txt")
    p.add_argument("--sigdir", default="gtdb95-evolpaths/signatures")
    p.add_argument("--sig-extension", default=".sig")
    p.add_argument("--alphabet", default="protein")
    p.add_argument("--ksize", default=10, type=int)
    p.add_argument("--scaled", default=100, type=int)
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
