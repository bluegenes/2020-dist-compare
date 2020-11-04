##
# Script to count and plot the number of shared k-mers at a taxonomic rank.
# modified from lca rankinfo.
##

import sys
import argparse
import csv
from collections import defaultdict
from collections import Counter

import sourmash
from sourmash.lca import lca_utils
#from sourmash.lca.command_rankinfo import make_lca_counts


def count_shared_kmers(dblist, min_num=2):
    """
    Collect counts of all the LCAs in the list of databases.
    CTB this could usefully be converted to a generator function.
    """

    # gather all hashvalue assignments from across all the databases
    assignments = defaultdict(set)
    for lca_db in dblist:
        for hashval, idx_list in lca_db.hashval_to_idx.items():
            if min_num and len(idx_list) < min_num:
                continue

            for idx in idx_list:
                lid = lca_db.idx_to_lid.get(idx)
                if lid is not None:
                    lineage = lca_db.lid_to_lineage[lid]
                    assignments[hashval].add(lineage)

    # find kmers shared across >1 lineage at each rank
    shared_kmers_by_lineage = defaultdict(int)
    for hashval, lineages in assignments.items():
        for rank in lca_utils.taxlist(strain=False).reverse():
             import pdb;pdb.set_trace()
             lin_at_rank = [lin.pop_to_rank(rank) for lin in lineages]
             shared_lineages = [lin for lin,count in Counter(lin_list).items() if count>=min_num]
             for lin in shared_lineages:
                 shared_kmers_by_lineage[lin] += 1

    return shared_kmers_by_lineage


def main(args):
    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.ksize, args.scaled)

    # count all the shared kmers across these databases
    counts = count_shared_kmers(dblist)

    # write out
    with open(args.csv, 'wt') as fp:
        hashes_by_lineage = csv.writer(fp)
        hashes_by_lineage.writerow(['rank', 'lineage', 'num_shared_kmers'])

        for lineage, shared_kmer_count in counts_by_lineage.items():
            rank = lineage[-1].rank
            hashes_by_lineage.writerow([rank, lca_utils.display_lineage(lineage), str(shared_kmer_count)])



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
