import os
import sys
import argparse
import glob
import pprint

import pandas as pd

from collections import defaultdict, namedtuple

anchorCompareM = namedtuple('anchorCompareM',
                           'comparison_name, anchor_name, ref_name, comparison_group, lowest_common_rank, mean_aai, std_aai, genes_in_anchor, genes_in_ref, orthologous_genes, orthologous_fraction')

compare_ranklist = ["genus", "family", "order", "class", "order", "class", "phylum", "superkingdom"]

def main(args):
    # get basename for these sequences
    compareInfo = pd.read_csv(args.comparison_info)
    compareInfo = compareInfo.groupby(["anchor_acc", "lowest_common_rank", "anchor_sciname"]).agg({"compare_accs":lambda x: list(x)})
    # load compareM files
    compareM_info= [tuple(x.strip().split(',')) for x in open(args.comparem_tsv_filecsv, "r")]
    anchor_results = []
    for (comparison_group, inF) in compareM_info:
        lowest_common_rank = comparison_group.split("-")[0]
        anchor_acc = comparison_group.split("-anchor")[1]

        aai_tsv = pd.read_csv(inF, sep = "\t", header=0)
        aai_tsv["#Genome A"] = aai_tsv["#Genome A"].str.rsplit("_", 1, expand=True)[0].str.replace('RS_' , '').str.replace('GB_' , '')
        aai_tsv["Genome B"] = aai_tsv["Genome B"].str.rsplit("_", 1, expand=True)[0].str.replace('RS_' , '').str.replace('GB_' , '')

        # first, let's get just the anchor results:
        # subset to just comparisons that include the anchor species ( later, we could pull out the all-by-all comparison if desired)
        anchor_compare_only = aai_tsv.loc[(aai_tsv["#Genome A"] == anchor_acc) | (aai_tsv["Genome B"] == anchor_acc)]

        # now loop through comparison accessions
        compare_accs = compareInfo.loc[(w.anchor, w.lcrank)]["compare_accs"].values[0]

        for compare_acc in compare_accs:
            comparison_name = f"{anchor_acc}_x_{compare_acc}"
            # pull out anchor x comparison (no guarantee on which column they'll be in)
            comparison_info = anchor_compare_only.loc[(aai_tsv["#Genome A"] == compare_acc) | (aai_tsv["Genome B"] == compare_acc)]
            if comparison_info["#Genome A"].values[0] == compare_acc:
                ref_genes = comparison_info["Genes in A"].values[0]
                anchor_genes = comparison_info["Genes in B"].values[0]
            else:
                ref_genes = comparison_info["Genes in B"].values[0]
                anchor_genes = comparison_info["Genes in A"].values[0]

            mean_aai = comparison_info["Mean AAI"].values[0]
            std_aai = comparison_info["Std AAI"].values[0]
            num_orthologous_genes = comparison_info["# orthologous genes"].values[0]
            orthologous_fraction = comparison_info["Orthologous fraction (OF)"].values[0]
            this_info = anchorCompareM(comparison_name, anchor_acc, compare_acc, comparison_group, lowest_common_rank, mean_aai, std_aai, anchor_genes, ref_genes, num_orthologous_genes, orthologous_fraction)
            anchor_results.append(this_info)

    # convert path aai comparison info to pandas dataframe
    anchor_aaiDF = pd.DataFrame.from_records(anchor_results, columns = anchorCompareM._fields)

    # print to csv
    anchor_aaiDF.to_csv(args.output_csv, index=False)
    print(f"done! path comparison info written to {args.output_csv}")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--comparem-tsv-filecsv")
    p.add_argument("--comparison-info")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
