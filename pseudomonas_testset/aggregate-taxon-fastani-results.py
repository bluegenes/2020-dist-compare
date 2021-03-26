import os
import sys
import argparse
import glob
import pprint

import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple

anchorfastANI = namedtuple('anchorFastANI',
                           'comparison_name, anchor_name, ref_name, comparison_group, lowest_common_rank, fastani_ident, num_bidirectional_fragment_mappings, total_query_fragments')

#compare_ranklist = ["genus", "family", "order", "class", "order", "class", "phylum", "superkingdom"]

def main(args):
    anchor_results = []

    compareInfo = pd.read_csv(args.comparison_info)
    compareInfo = compareInfo.groupby(["anchor_acc", "lowest_common_rank", "anchor_sciname"]).agg({"compare_accs":lambda x: list(x)})
    # loop through fastani comparison files
    fastani_info= [tuple(x.strip().split(',')) for x in open(args.fastani_filecsv, "r")]
    for (comparison_group, inF) in fastani_info:
        lowest_common_rank = comparison_group.split("-")[0]
        anchor_acc = comparison_group.split("-anchor")[1]

        fastani = pd.read_csv(inF, sep = "\t", header=None, names=['anchor','ref','fastani_ident','count_bidirectional_frag_mappings','total_query_frags'])
        if fastani.empty:
            continue
        fastani["anchor"] = fastani["anchor"].str.rsplit("/", 1, expand=True)[1].str.rsplit("_genomic.fna.gz", 1, expand=True)[0]
        fastani["ref"] = fastani["ref"].str.rsplit("/", 1, expand=True)[1].str.rsplit("_genomic.fna.gz", 1, expand=True)[0]
        fastani.set_index("ref",inplace=True)

        # now check comparisons for results:
        compare_accs = compareInfo.loc[(anchor_acc, lowest_common_rank)]["compare_accs"].values[0]

        for compare_acc in compare_accs:
            comparison_name = f"{anchor_acc}_x_{compare_acc}"
            fastani_ident, bidirectional_mappings, query_frags = np.nan, np.nan, np.nan
            if compare_acc in fastani.index:
                fastani_ident = fastani.at[compare_acc, "fastani_ident"]
                bidirectional_mappings = fastani.at[compare_acc, "count_bidirectional_frag_mappings"]
                query_frags = fastani.at[compare_acc, "total_query_frags"]

            this_info = anchorfastANI(comparison_name, anchor_acc, compare_acc, comparison_group, lowest_common_rank, fastani_ident, bidirectional_mappings, query_frags)
            anchor_results.append(this_info)

    # convert path aai comparison info to pandas dataframe
    anchor_aniDF = pd.DataFrame.from_records(anchor_results, columns = anchorfastANI._fields)

    # print to csv
    anchor_aniDF.to_csv(args.output_csv, index=False)
    print(f"done! path comparison info written to {args.output_csv}")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--fastani-filecsv")
    p.add_argument("--comparison-info")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
