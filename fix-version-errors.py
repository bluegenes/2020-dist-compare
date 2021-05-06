import os
import sys
import argparse
import glob
import pprint

import pandas as pd

def main(args):
    metadataDF = pd.read_csv(args.gtdb_metadata, header=0, low_memory=False)
    metadataDF["accession"] = metadataDF["accession"].str.replace("GB_", "").str.replace("RS_", "")
    # read in info about version errors encountered during downloads
    version_err = pd.read_csv(args.version_errors_log, header=None, names = ["accession", "prev_version", "full_name"])
    verrD = pd.Series(version_err["prev_version"].values,index=version_err["accession"]).to_dict()
   # replace accessions to new accs
    metadataDF["accession"].replace(verrD, inplace=True)
    import pdb;pdb.set_trace()
    metadataDF.to_csv(args.output_csv, index=False)


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--gtdb-metadata",  default = "/group/ctbrowngrp/gtdb/gtdb-r202.metadata.csv.gz")
    p.add_argument("--version-errors-log", default="/home/ntpierce/sourmash_databases/output.gtdb-databases/accession_version_errors.csv")
    p.add_argument("--output-csv", default="/group/ctbrowngrp/gtdb/gtdb-r202.metadata.v2.csv.gz")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

