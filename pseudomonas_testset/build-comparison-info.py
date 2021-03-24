import os
import sys
import argparse
import glob
import pprint

import pandas as pd

from collections import defaultdict, namedtuple


# build long form comparison csv, each line containing the anchor acc and the comparison genome.
# use these later to build appropriate filelists for fastani, compareM and sourmash comparisons
# AND/OR ... build and dump a json!? easier!

ComparisonInfo = namedtuple('ComparisonInfo',
                          'anchor_acc, lowest_common_rank, anchor_sciname, tax_id, compare_accs')

def main(args):
    infoDF = pd.read_csv(args.taxinfo_csv)
    comparison_ranks = ["genus", "species", "strain"]
    cInfo = []
    #cJ = {}
    for c_rank in comparison_ranks:
        # genus comparisons
        if c_rank == "genus":
            # pick first with Pseudomonas parent_tax_id
            psd_taxid = 286
            genus_anchor_acc = infoDF[infoDF["parent_tax_id"] == 286]["accession"][0]
            genus_name = infoDF[infoDF["accession"] == genus_anchor_acc]["sci_name"].tolist()[0]
            genus_compare_accessions = infoDF["accession"].tolist()
            genus_compare_accessions.remove(genus_anchor_acc)
            #cInfo.append(ComparisonInfo("genus", genus_name, psd_taxid, genus_anchor_acc, genus_compare_accessions))
            #cJ[genus_anchor_acc] = ("genus", genus_name, psd_taxid, genus_compare_accessions)
            for gc_acc in genus_compare_accessions:
                #multiple rows per comparison (long format)
                cInfo.append(ComparisonInfo(genus_anchor_acc, "genus", genus_name, psd_taxid, gc_acc))
        # species comparisons
        if c_rank == "species":
            # for species, don't pick 0th as anchor, pick 1th. Bc can only pairwise compare things with >1 genome
            species_anchors=infoDF.groupby(["tax_id"], as_index=False).nth(1)[["accession", "tax_id"]].to_records(index=False).tolist()
            for (sp_anchor_acc, taxid) in species_anchors:
                species_name = infoDF[infoDF["accession"] == sp_anchor_acc]["sci_name"].tolist()[0]
                species_acc = infoDF[infoDF["tax_id"] == taxid]["accession"].tolist()
                species_acc.remove(sp_anchor_acc)
                #cInfo.append(ComparisonInfo("species", species_name, taxid, sp_anchor_acc, species_acc))
                for sp_acc in species_acc:
                    cInfo.append(ComparisonInfo(sp_anchor_acc, "species", species_name, taxid, sp_acc))

        # strain comparisons
        if c_rank == "strain":
            strain_anchors=infoDF.groupby(["strain"], as_index=False).nth(1)[["accession", "strain"]].to_records(index=False).tolist()
            for (strain_anchor_acc, strain) in strain_anchors:
                strain_name = infoDF[infoDF["accession"] == strain_anchor_acc]["sci_name"].tolist()[0]
                strain_acc = infoDF[infoDF["strain"] == strain]["accession"].tolist()
                strain_acc.remove(strain_anchor_acc)
                #cInfo.append(ComparisonInfo("strain", strain_name, strain, strain_anchor_acc, strain_acc))
                for st_acc in strain_acc:
                    cInfo.append(ComparisonInfo(strain_anchor_acc, "strain", strain_name, strain, st_acc))


    # convert to pandas DF and write csv:
    compareDF = pd.DataFrame.from_records(cInfo, columns = ComparisonInfo._fields)
    # sort by comparison level
    compareDF.sort_values("lowest_common_rank", inplace=True)
    compareDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--taxinfo-csv",  default = "pseudomonas_genomes.info.CompleteGenome.lineages.csv")
    p.add_argument("--output-csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

