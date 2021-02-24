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
                           'comparison_name, comparison_level, alphabet, ksize, scaled, jaccard, max_containment, anchor_hashes, query_hashes, num_common')

def compare_sigs(sigA, sigB, comparison_name, comparison_level, alpha, ksize, scaled):
    sigA_numhashes = len(sigA.minhash.hashes)
    sigB_numhashes = len(sigB.minhash.hashes)
    intersect_numhashes = sigA.minhash.count_common(sigB.minhash)
    jaccard = sigA.jaccard(sigB)
    contain1 = sigA.contained_by(sigB)
    contain2 = sigB.contained_by(sigA)
    max_contain = max(contain1,contain2)
    return CompareResult(comparison_name, comparison_level, alpha, ksize, scaled, jaccard, max_contain, sigA_numhashes, sigB_numhashes, intersect_numhashes)

def main(args):
    sigdir = args.sigdir
    alpha = args.compare_alphabet
    if alpha == "nucleotide":
        moltype = "DNA"
    else:
        moltype = alpha
    print(moltype)
    ksize = args.ksize
    scaled = args.scaled
    ext = args.sig_ext

    infoDF = pd.read_csv(args.taxinfo_csv)
    # just do species-level right now. Some strains have multiple genomes per strain --> HANDLE LATER
    #all_strains = taxinfoDF["strain"].unique()
    #strainDF = queryDF.copy().groupby([rep_col], as_index=False).nth(select_n)
    #strainDF.dropna(inplace=True)
    ## tax_id == species. parent_tax_id !=genus though... just pick an anchor and compare all against it
    # for species, don't pick 0th as anchor, pick 1th. Bc can only pairwise compare things with >1 genome
    species_anchors =infoDF.groupby(["tax_id"], as_index=False).nth(1)[["accession", "tax_id"]].to_records(index=False).tolist()

    genus_anchor = infoDF["accession"][0]
    # get genus anchor sigfile
    genus_anchorSigF = os.path.join(sigdir, genus_anchor + ext)
    # all accessions to compare to genus sig
    all_accessions = infoDF["accession"].tolist()
    all_accessions.remove(genus_anchor)

    compareInfo=[]
    ### compare sigs at each ksize
    print(f"Working on alphabet: {alpha}, ksize: {ksize} (selecting ksize: {ksize})")

    # load genus anchor
    genus_sig_gen = sourmash.sourmash_args.load_file_as_signatures(genus_anchorSigF,ksize=ksize, select_moltype=moltype)
    genus_sig = next(genus_sig_gen)

    for (sp_anchor, taxid) in species_anchors:
        species_name = infoDF[infoDF["accession"] == sp_anchor]["sci_name"].tolist()[0]
        species_acc = infoDF[infoDF["tax_id"] == taxid]["accession"].tolist()
        print(f"Working on taxid {taxid}, sci name: {species_name} ({len(species_acc)} genomes)")
        species_acc.remove(sp_anchor)
        # get species anchor sigfile
        anchor_sigF = os.path.join(sigdir, sp_anchor + ext)
        # get list of accessions to compare
        # load species anchor
        #species_anchor_sig = list(sourmash.signature.load_signatures(anchor_sigF, ksize=ksize, select_moltype=moltype))[0]
        species_anchor_gen = sourmash.sourmash_args.load_file_as_signatures(anchor_sigF,ksize=ksize, select_moltype=moltype)
        species_anchor_sig = next(species_anchor_gen)

        # compare
        for acc in species_acc:
            full_sigF = os.path.join(sigdir, acc + ext)
            query_gen = sourmash.sourmash_args.load_file_as_signatures(full_sigF, ksize=ksize, select_moltype=moltype)
            query_sig = next(query_gen)
            #query_sig = list(sourmash.signature.load_signatures(full_sigF, ksize=ksize, select_moltype=moltype))[0]
            # first, do species-level comparison
            comparison_name = sp_anchor + "__x__" + acc
            sComparison = compare_sigs(species_anchor_sig, query_sig, comparison_name, "species", alpha, ksize, scaled)
            compareInfo.append(sComparison)
            # now, compare against genus anchor sig
            comparison_name = genus_anchor + "__x__" + acc
            gComparison = compare_sigs(genus_sig, query_sig, comparison_name, "genus", alpha, ksize, scaled)
            compareInfo.append(gComparison)
            # now this comparison is done, remove from all_accessions
            all_accessions.remove(acc)

    # now do remaining genus comparisons at this ksize/scaled
    print(f"Now running genus-level comparisons for remaining {len(all_accessions)} genomes")
    for acc in all_accessions:
        full_sigF = os.path.join(sigdir, acc + ext)
        query_gen = sourmash.sourmash_args.load_file_as_signatures(full_sigF, ksize=ksize, select_moltype=moltype)
        query_sig = next(query_gen)
        #query_sig = list(sourmash.signature.load_signatures(full_sigF, ksize=ksize, select_moltype=moltype))[0]
        # now, compare against genus anchor sig
        comparison_name = genus_anchor + "__x__" + acc
        gComparison = compare_sigs(genus_sig, query_sig, comparison_name, "genus", alpha, ksize, scaled)
        compareInfo.append(gComparison)

    # convert to pandas DF and write csv:
    compareDF = pd.DataFrame.from_records(compareInfo, columns = CompareResult._fields)
    # sort by comparison level
    compareDF.sort_values("comparison_level", inplace=True)
    compareDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--taxinfo-csv",  default = "pseudomonas_genomes.info.CompleteGenome.lineages.csv")
    p.add_argument("--sigdir", default = "")
    p.add_argument("--sig-ext", default = ".genomic.sig")
    p.add_argument("--compare-alphabet", default="nucleotide")
    p.add_argument("--scaled", default=31, type=int)
    p.add_argument("--ksize", default=1000, type=int)
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

