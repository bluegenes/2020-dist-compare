##
# Script to count and plot the number of shared k-mers at a taxonomic rank.
# modified from lca rankinfo.
##

import sys
import argparse
import csv
from collections import defaultdict

import sourmash
from sourmash.lca import lca_utils
from sourmash.lca.command_rankinfo import make_lca_counts

def main(args):
    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.ksize, args.scaled)

    # count all the LCAs across these databases
    counts = make_lca_counts(dblist)

    # collect counts across all ranks
    counts_by_rank = defaultdict(int)
    counts_by_lineage = defaultdict(int)
    for lineage, count in counts.items():
        if lineage:
            lineage_tup = lineage[-1]
            counts_by_rank[lineage_tup.rank] += count
            # count by lineage instead of just rank
            counts_by_lineage[lineage] += count

    # output!
    total = float(sum(counts_by_rank.values()))
    if total == 0:
        notify("(no hashvals with lineages found)")
    else:
        for rank in lca_utils.taxlist():
            count = counts_by_rank.get(rank, 0)
            print('{}: {} ({:.1f}%)'.format(rank, count, count / total * 100.))


        # counts_by_lineage now has lca kmer count. BUT, shared k-mer count = lca kmer count + lca kmer count at all lower levels (e.g. for a genus underneath an order).
        # now output counts_by_lineage to csv: lineage, shared (lca) kmer count
        with open(args.csv, 'wt') as fp:
            hashes_by_lineage = csv.writer(fp)
            hashes_by_lineage.writerow(['rank', 'lineage', 'num_lca_kmers'])
            # Could sort lineages by num lca kmers... meh.
            #counts_by_lineage_items = list(counts_by_lineage.items())
            #counts_by_lineage_items.sort(...)
            for lineage, lca_count in counts_by_lineage.items():
                rank = lineage[-1].rank
                hashes_by_lineage.writerow([rank, lca_utils.display_lineage(lineage), str(lca_count)])



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, nargs="+")
    p.add_argument("--ksize", type=int)
    p.add_argument("--scaled")
    p.add_argument("--csv")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)
