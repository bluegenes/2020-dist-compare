import os
import sys
import argparse
import glob
import pprint

import pandas as pd
from sourmash import lca_utils
from collections import namedtuple

ComparisonInfo = namedtuple('ComparisonInfo',
                          'anchor_acc, lowest_common_rank, lowest_common_taxon, compare_accs')



def select_representative_genomes(rDF, ranklist = lca_utils.taxlist(include_strain=False)):
    anchors_at_ranks = {}
    for rank in ranklist:
        reps = rDF.groupby(rank, as_index=False).nth(0).index
        anchors_at_ranks[rank] = reps
    return anchors_at_ranks

def taxonomy_at_ranks(row):
    rank_lineages = []
    lin = lca_utils.make_lineage(row[taxonomy])
    for rank in lca_utils.taxlist(include_strain=False):
        rank_lin = lca_utils.pop_to_rank(lin, rank)
        rank_lin = ";".join(lca_utils.zip_lineage(rank_lin, include_strain=False, truncate_empty=True))
        row[rank] = rank_lin
    return row


def main(args):
    metadataDF = pd.read_csv(args.gtdb_metadata, header=0, low_memory=False)
    metadataDF["accession"] = metadataDF["accession"].str.replace("GB_", "").str.replace("RS_", "")
    taxonomy = "gtdb_taxonomy"
    # subset to just gtdb relevant taxonomy information
    if args.taxonomy == "gtdb":
        mDF = metadataDF[["accession", "gtdb_taxonomy"]].copy()
        mDF = mDF.apply(taxonomy_at_ranks,axis=1)
    elif args.taxonomy == "ncbi":
        mDF = metadataDF[["accession", "ncbi_taxonomy"]].copy()
        mDF = mDF.apply(taxonomy_at_ranks,axis=1)
        taxonomy = "ncbi_taxonomy"
    mDF["signame"] = mDF["accession"] + " " + mDF["species"].str.replace("s__", "")
    mDF.set_index("accession", inplace=True)

     # select anchors!
    if args.taxonomy == "gtdb" and not args.randomly_select_anchors:
       # use to the representative genomes, only relevant for gtdb?
        representatives_only = metadataDF[metadataDF["gtdb_representative"] == "t"]
        reps_only = mDF[(mDF.index).isin(representatives_only["accession"])]
        anchors_at_ranks = select_representative_genomes(reps_only)
    else:
        # select randomly
        anchors_at_ranks = select_representative_genomes(mDF)

    # for each anchor, find the accessions we need to compare to
    cInfo = []
    for i, rank in enumerate(lca_utils.taxlist(include_strain=False)):
        # 0 - superk, 1-phylum, etc
        for acc in anchors_at_ranks[rank]:
            full_taxonomy = mDF.at[acc, taxonomy]

            taxon_accessions = list(mDF[mDF[rank] == common_taxon].index)

            # do I know for sure this is the lowest common taxon?? What situations might it _not_ be?

            ### working here -- now that I'm using full taxonomy to the rank of interest, these should be right?? CHECK.
            fo acc in taxon_accessons:
                if mDF.at[acc, prev_rank] !=
            acc_to_compare = taxon_accessions
            # if there's only one genome in this group, nothing to compare!
            if len(acc_to_compare) == 1:
                continue
             # remove anchor accession
            acc_to_compare.remove(acc)
            acc_to_compare = ";".join(acc_to_compare)
            cInfo.append(ComparisonInfo(acc, rank, lowest_common_taxon, acc_to_compare))

    # convert to pandas DF and write csv:
    compareDF = pd.DataFrame.from_records(cInfo, columns = ComparisonInfo._fields)
    # sort by comparison level
    compareDF.sort_values("lowest_common_rank", inplace=True)
    compareDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--gtdb-metadata",  default = "/group/ctbrowngrp/gtdb/gtdb-r202.metadata.csv.gz")
    p.add_argument("--taxonomy", default="gtdb", choices = ["gtdb", "ncbi"])
    p.add_argument("--randomly-select-anchors", action="store_true")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

