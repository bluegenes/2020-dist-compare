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

def pop_to_rank(lin, rank):
    "Remove lineage tuples from given lineage `lin` until `rank` is reached."
    lin = list(lin)

    txl = lca_utils.taxlist()
    before_rank = []
    for txl_rank in txl:
        if txl_rank != rank:
            before_rank.append(txl_rank)
        else:
            break

    # are we already above rank?
    if lin and lin[-1].rank in before_rank:
        return tuple(lin)

    while lin and lin[-1].rank != rank:
        lin.pop()

    return tuple(lin)



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
    rev_taxlist = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']
    shared_kmers_by_lineage = defaultdict(int)
    total_kmers = len(assignments.keys())
    for hashval, lineages in assignments.items():
        for rank in rev_taxlist:
            lin_at_rank = [pop_to_rank(lin, rank) for lin in lineages]
            shared_lineages = [lin for lin,count in Counter(lin_at_rank).items() if count>=min_num]
            import pdb;pdb.set_trace()
            for lin in shared_lineages:
                shared_kmers_by_lineage[lin] += 1

    return shared_kmers_by_lineage, total_kmers


def main(args):
    # load all the databases
    dblist, ksize, scaled = lca_utils.load_databases(args.db, args.ksize, args.scaled)

    # count all the shared kmers across these databases
    counts, total_kmer_count= count_shared_kmers(dblist)

    # write out
    with open(args.csv, 'wt') as fp:
        hashes_by_lineage = csv.writer(fp)
        hashes_by_lineage.writerow(['rank', 'lineage', 'num_shared_kmers', 'percent_shared_kmers'])

        for lineage, shared_kmer_count in counts.items():
            rank = lineage[-1].rank
            percent_shared_kmers = float(shared_kmer_count)/total_kmer_count
            hashes_by_lineage.writerow([rank, lca_utils.display_lineage(lineage), str(shared_kmer_count), str(round(percent_shared_kmers, 2))])



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
