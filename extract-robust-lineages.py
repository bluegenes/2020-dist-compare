#! /usr/bin/env python3
import argparse
import sourmash.lca
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca import LineagePair
import sys
import collections
import csv
import random


def pop_to_rank(lin, rank):
    "Remove lineage tuples from given lineage `lin` until `rank` is reached."
    lin = list(lin)

    txl = sourmash.lca.lca_utils.taxlist()
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


def is_lineage_match(lin_a, lin_b, rank):
    """
    check to see if two lineages are a match down to given rank.
    """
    for a, b in zip(lin_a, lin_b):
        assert a.rank == b.rank
        if a.rank == rank:
            if a == b:
                return 1
        if a != b:
            return 0

    return 0


def get_name_at_rank(lineage, rank):
    for pair in lineage:
        if pair.rank == rank:
            return pair.name

    raise Exception


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lineage_csv')
    p.add_argument('output_filenames')
    p.add_argument("pathinfo")
    p.add_argument('-N', '--num_paths', type=int, default=0)
    args = p.parse_args()


    assignments, n_rows = load_taxonomy_assignments(args.lineage_csv,
                                                    start_column=4)
    print(f'loaded {len(assignments)} assignments.')

    # here, 'assignments' is a dictionary where the keys are the
    # identifiers in column 1, and the values are tuples of LineagePairs.

    paths_to_idents = collections.defaultdict(set)

    with open(args.pathinfo, "w") as outP:
        header=["accession", "path", "rank", "lineage"]
        outP.write("\t".join(header) + "\n")
        # connect every lineage in lineage_paths to their identifiers.
        for ident, lineage in assignments.items():
            # bring back to each rank, skipping species
            for rank in sourmash.lca.taxlist(include_strain=False):
                if rank == 'species': continue

                tup = pop_to_rank(lineage, rank)
                assert tup
                assert tup[-1].rank == rank

                # record identifer for that lineage.
                paths_to_idents[tup].add(ident)

        # now, collect all genus level paths
        genus_paths = set()
        for k, v in paths_to_idents.items():
            if k[-1].rank == 'genus':
                genus_paths.add(k)

        # for each genus level path, pull out actual accessions for the
        # entire path.
        genus_paths = list(sorted(genus_paths))
        extract_idents = set()

        paths_found = 0
        for chosen_genus in genus_paths:
            if args.num_paths and paths_found == args.num_paths:
                break

            # for this chosen genus, find identifiers for lineages that differ at
            # each step. lineages may not have this, note - CTB explain later :)
            d = {}
            lineage = list(chosen_genus)
            assert lineage[-1].rank == 'genus'
            last_rank = 'species'

            # pick a specific identifier and grab full lineage for it:
            chosen_ident = next(iter(paths_to_idents[chosen_genus]))
            chosen_lineage = assignments[chosen_ident]

            # now find differences at every level.
            n_found = 0
            while lineage:
                rank = lineage[-1].rank
                idents = paths_to_idents[tuple(lineage)]

                for ident in idents:
                    this_lineage = assignments[ident]

                    # find ones that match at this level, and not previous level
                    if is_lineage_match(this_lineage, chosen_lineage, rank) and \
                       not is_lineage_match(this_lineage, chosen_lineage, last_rank):
                         d[rank] = ident
                         n_found += 1
                         break

                last_rank = rank
                lineage.pop()

            if n_found < 6:
                continue                      # skip on to next genus-level path

            assert n_found == 6

            paths_found += 1
            pathname = f"path{str(paths_found)}" # use path number as a unique identifier for each path
            ### now, start tracking down identifiers to extract!

            extract_idents.add(chosen_ident)

            # print some stuff out --
            print('------------- path:', paths_found)
            for rank in sourmash.lca.taxlist():
                if rank == 'species':
                    print(' ', "species", chosen_ident, sourmash.lca.display_lineage(chosen_lineage))
                    outP.write("\t".join([chosen_ident, pathname, "species", sourmash.lca.display_lineage(chosen_lineage)]) + "\n")
                    break

                assert rank in d
                ident = d[rank]
                this_lineage = assignments[ident]
                print(' ', rank, ident, sourmash.lca.display_lineage(this_lineage))
                outP.write("\t".join([ident, pathname, rank, sourmash.lca.display_lineage(this_lineage)]) + "\n")

                extract_idents.add(ident)

        # last but not least -- extract filenames
        with open(args.output_filenames, 'wt') as fp:
            r = csv.DictReader(open(args.lineage_csv, 'rt'))
            for row in r:
                if row['accession'] in extract_idents:
                    filename = row['filename']
                    print(filename, file=fp)


if __name__ == '__main__':
    main()
