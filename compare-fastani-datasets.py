import os
import sys
import argparse
import glob
import pprint

import sourmash
import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple

#def convert_to_pdist(corr, ksize):
#    # return proportion of observed differences
#    if corr ==0:
#        return 1.0 # what should this be? (jaccard of 1.0 returns 0)
#    p = 1 - np.power(2*corr/(corr + 1),(1/float(ksize)))
#    return p
CompareResult = namedtuple('CompareResult',
                           'signame, moltype, ksize, scaled, jaccard, max_containment, anchor_hashes, query_hashes, num_common')


def main(args):
    # get basename for these sequences
    if args.input_alphabet == "nucleotide":
        alphabets_to_compare = ["nucleotide"] #, "protein", "dayhoff", "hp"]
    else:
        alphabets_to_compare = ["protein", "dayhoff", "hp"]
    ksize_info = {"nucleotide": [21,31,51], "protein": [7,8,9,10,11,12], "dayhoff": [15,16,17,18,19], "hp": [33,35,37,39,42]}
    scaled_info = {"nucleotide": [1000], "protein": [100], "dayhoff": [100], "hp": [100]}

    #compareInfo=defaultdict(list)
    compareInfo=[]
    # get sigfile names and anchor sigfile
    sigfiles = [x.strip() for x in open(args.siglist, "r")]
    anchor_sigF = args.anchor_sig
    # remove anchorsig from info_filenames
    sigfiles.remove(anchor_sigF)
    anchor_sigfile = os.path.join(args.sig_dir, anchor_sigF)

    ### compare sigs at each ksize
    translate = False
    for alpha in alphabets_to_compare:
        if alpha == "nucleotide":
            moltype = "dna"
        else:
            moltype = alpha
        for ksize in ksize_info[alpha]:
            select_ksize = ksize
            if moltype != "dna":
                select_ksize = int(ksize)*3
            print(f"Working on alphabet: {alpha}, ksize: {ksize} (selecting ksize: {select_ksize})")
            for scaled in scaled_info[alpha]:
                anchor_sig = list(sourmash.signature.load_signatures(anchor_sigfile, ksize=select_ksize, select_moltype=moltype))[0]
                anchor_hashes = len(anchor_sig.minhash.hashes)
                for sigF in sigfiles:
                    full_sigF = os.path.join(args.sig_dir, sigF)
                    query_sig = list(sourmash.signature.load_signatures(full_sigF, ksize=select_ksize, select_moltype=moltype))[0]
                    query_hashes = len(query_sig.minhash.hashes)
                    intersect_hashes = anchor_sig.minhash.count_common(query_sig.minhash)
                    jaccard = anchor_sig.jaccard(query_sig)
                    contain1 = anchor_sig.contained_by(query_sig)
                    contain2 = query_sig.contained_by(anchor_sig)
                    max_contain = max(contain1,contain2)
                    compareInfo.append(CompareResult(os.path.basename(sigF), alpha, ksize, scaled, jaccard, max_contain, anchor_hashes, query_hashes, intersect_hashes))

    # convert to pandas DF and write csv:
    infoDF = pd.DataFrame.from_records(compareInfo, columns = CompareResult._fields)
    infoDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--siglist")
    p.add_argument("--sig-dir", default = "")
    p.add_argument("--anchor-sig")
    p.add_argument("--input-alphabet", default="nucleotide")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

