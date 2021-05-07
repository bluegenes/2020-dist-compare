import os
import sys
import argparse
import glob
import pprint
import pytest

import pandas as pd
from sourmash.lca import lca_utils
from collections import namedtuple
from itertools import combinations

ComparisonInfo = namedtuple('ComparisonInfo',
                          'comparison_name, accessionA, accessionB, lowest_common_rank, lowest_common_lineage, is_representativeA, is_representativeB')


def find_lca(linA,linB, reverse_taxlist=["species", "genus", "family", "order", "class", "phylum", "superkingdom"]):
    lca_rank=None
    lca_lin=None
    for rank in reverse_taxlist:
            # do we have a lca lineage match at this rank?
        if lca_utils.is_lineage_match(linA, linB, rank):
            lca_rank = rank
            lca_lin = lca_utils.pop_to_rank(linA, lca_rank)
            break
    return lca_rank, lca_lin



def main(args):
    metadataDF = pd.read_csv(args.gtdb_metadata, header=0, low_memory=False)
    metadataDF["accession"] = metadataDF["accession"].str.replace("GB_", "").str.replace("RS_", "")
    metadataDF.set_index("accession", inplace=True)
    representative_accs = list(metadataDF[metadataDF["gtdb_representative"] == "t"].index) # does this need to be a list?
    # which taxonomy?
    taxonomy = "gtdb_taxonomy"
    if args.taxonomy == "ncbi":
        taxonomy = "ncbi_taxonomy"

    cInfo = []
    # get dictionary of species_path :: accession
    acc2path = pd.Series(metadataDF[taxonomy].values,index=metadataDF.index).to_dict()
    reverse_taxlist = list(lca_utils.taxlist(include_strain=False))[::-1]
    # iterate through all possible pairwise comparisons
    for n, (acc1, acc2) in enumerate(combinations(acc2path.keys(), r = 2)):
        comparison_name = f"{acc1}_x_{acc2}"
        if n !=0 and n % 10000 == 0:
            print(f"... checking {n}th comparison, {comparison_name}\n")
        # check if either are gtdb representative genomes
        rep1,rep2=False,False
        if acc1 in representative_accs:
          rep1 = True
        if acc2 in representative_accs:
          rep2 = True
        lin1 = lca_utils.make_lineage(acc2path[acc1])
        lin2 = lca_utils.make_lineage(acc2path[acc2])
        lca_rank,lca_lin = find_lca(lin1,lin2)
       # if there's an LCA rank, store the comparison
        if lca_rank:
            lca_lin = ";".join(lca_utils.zip_lineage(lca_lin, truncate_empty=True))
            if n !=0 and n % 1000 == 0:
                print(f"...adding lca comparison for rank {lca_rank}: {comparison_name}")
            cInfo.append(ComparisonInfo(comparison_name,acc1,acc2,lca_rank,lca_lin,rep1,rep2))

    # convert to pandas DF and write csv:
    compareDF = pd.DataFrame.from_records(cInfo, columns = ComparisonInfo._fields)
    # sort by comparison level
    compareDF.sort_values("lowest_common_rank", inplace=True)
    compareDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--gtdb-metadata",  default = "/group/ctbrowngrp/gtdb/gtdb-r202.metadata.v2.csv.gz")
    p.add_argument("--taxonomy", default="gtdb", choices = ["gtdb", "ncbi"])
    p.add_argument("--output-csv", default="gtdb-pairwise-lca-comparisons.csv.gz")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)


def test_build_comparison_lowest_lca():
    lin1 = lca_utils.make_lineage('a;b;c')
    print(lin1)
    lin2 = lca_utils.make_lineage('a;b;c')

    lca_rank, lca_lin = find_lca(lin1,lin2)
    print("LCA RANK", lca_rank)
    print("LCA LIN", lca_lin)
    assert lca_rank == "class"
    assert lca_lin == lin1

def test_build_comparison_mid_lca():
    lin1 = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;b;d')
    true_lca_lin = lca_utils.make_lineage('a;b')

    lca_rank, lca_lin = find_lca(lin1,lin2)
    assert lca_rank == "phylum"
    assert lca_lin == true_lca_lin

def test_build_comparison_highest_lca():
    lin1 = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('a;d;e')
    true_lca_lin = lca_utils.make_lineage('a')

    lca_rank, lca_lin = find_lca(lin1,lin2)
    assert lca_rank == "superkingdom"
    assert lca_lin == true_lca_lin

def test_build_comparison_no_lca():
    lin1 = lca_utils.make_lineage('a;b;c')
    lin2 = lca_utils.make_lineage('d;e;f')

    lca_rank, lca_lin = find_lca(lin1,lin2)
    assert lca_rank == None
    assert lca_lin == None

def test_reverse_taxlist():
    reverse_taxlist = list(lca_utils.taxlist(include_strain=False))[::-1]
    assert reverse_taxlist == ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

