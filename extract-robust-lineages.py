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
    args = p.parse_args()

    assignments, n_rows = load_taxonomy_assignments(args.lineage_csv,
                                                    start_column=3)
    print(f'loaded {len(assignments)} assignments.')

    # here, 'assignments' is a dictionary where the keys are the
    # identifiers in column 1, and the values are tuples of LineagePairs.

    # let's start by filtering out all lineages where there is only one
    # at that rank/name.

    # build a tree. each node is a dictionary; keys are LineagePair,
    # values are dictionaries of next taxonomic level beneath.
    tree = sourmash.lca.build_tree(assignments.values())

    # recurse through and identify all lineages with at least two
    # children, all the way down to species. This should ensure
    # that there are at least two separate descendants of each node.
    filtered_lineages = set()

    def recurseme(node, path):
        # if this node has two or more children, keep on going.
        if len(node) >= 2:
            for k, v in node.items():
                downpath = list(path) + [k]   # track path => next level
                recurseme(v, downpath)
        elif path[-1].rank == 'species':      # whups, we are at end!
            # add to our set of interesting lineages
            filtered_lineages.add(tuple(path))
        else:                                 # ends not in species, or < 2
            pass

    recurseme(tree, [])

    # build another tree of just the filtered lineages --
    ftree = sourmash.lca.build_tree(filtered_lineages)

    # track all of the paths in this tree in a set.
    lineage_paths = set()
    def recurseme2(node, path):
        # for each node, record this path somewhere
        lineage_paths.add(tuple(path))

        # recurse!
        for k, v in node.items():
            downpath = list(path) + [k]   # track path => next level
            recurseme2(v, downpath)

    recurseme2(ftree, [])

    # next, build a dictionary of lineage_paths to identifiers
    # skipping species-level ranks.
    paths_to_idents = {}
    for path in lineage_paths:
        if path and path[-1].rank != 'species':
            paths_to_idents[path] = set()

    # connect every lineage in lineage_paths to their identifiers.
    for ident, lineage in assignments.items():
        # for every lineage that's in filtered lineages,
        if lineage in filtered_lineages:

            # bring back to each rank, skipping species
            for rank in sourmash.lca.taxlist(include_strain=False):
                if rank == 'species': continue

                tup = pop_to_rank(lineage, rank)
                assert tup
                assert tup in lineage_paths
                assert tup[-1].rank == rank

                # record identifer for that lineage.
                paths_to_idents[tup].add(ident)

    # double-check: should have two identifiers for each path!
    n_genus = 0
    for k, v in paths_to_idents.items():
        assert len(v) >= 2
        if k[-1].rank == 'genus':
            n_genus += 1

    print(f'found {n_genus} genus level pairs')

    # now, collect all genus level paths
    genus_paths = set()
    for k, v in paths_to_idents.items():
        if k[-1].rank == 'genus':
            genus_paths.add(k)

    # pick N genus-level paths --
    genus_paths = list(sorted(genus_paths))
    extract_idents = set()

    # for each genus level path, pull out some actual accessions for the
    # entire path.
    for chosen_genus in genus_paths:

        # for this chosen genus, find identifiers for lineages that differ at
        # each step. lineages may not have this, note - CTB explain later :)
        d = {}
        lineage = list(chosen_genus)
        assert lineage[-1].rank == 'genus'
        last_rank = 'species'

        # pick a specific identifier and grab full lineage for it:
        chosen_ident = next(iter(paths_to_idents[chosen_genus]))
        chosen_lineage = assignments[chosen_ident]
        extract_idents.add(chosen_ident)

        # now find differences at every level.
        while lineage:
            rank = lineage[-1].rank
            idents = paths_to_idents[tuple(lineage)]
            found = False

            for ident in idents:
                this_lineage = assignments[ident]
                
                # find ones that match at this level, and not previous level
                if is_lineage_match(this_lineage, chosen_lineage, rank) and \
                   not is_lineage_match(this_lineage, chosen_lineage, last_rank):
                     d[rank] = ident
                     found = True
                     break

            last_rank = rank
            lineage.pop()

        found = True
        for rank in sourmash.lca.taxlist():
            if rank == 'species':
                break

            if rank not in d:
                found = False

        if not found:
            continue                      # skip out
        

        print('-------------')
        for rank in sourmash.lca.taxlist():
            if rank == 'species':
                break

            assert rank in d
            ident = d[rank]
            this_lineage = assignments[ident]
            print(rank, ident)
            print(sourmash.lca.display_lineage(this_lineage))
            print('')

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
