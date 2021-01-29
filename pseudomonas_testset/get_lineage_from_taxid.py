import os
import sys
import glob
import argparse
import pandas as pd
import taxonomy
import sourmash

def get_row_taxpath(row, taxo, ranks, taxid_col="ncbi-taxid"):
    # modified from https://github.com/dib-lab/2019-12-12-sourmash_viz/blob/1223b736add63ea49108eecceb3f4bca85c78492/src/gather_to_opal.py
    current_taxid = str(row[taxid_col])
    try:
        current_rank = taxo.rank(current_taxid)
    except Exception as e:
        # pyo3 doesn't export Exceptions properly, so doing this for now...
        if "TaxonomyError" in repr(e):
            print(e)
            # a TaxonomyError is raised when the taxid is not in the taxonomy.
            # Skipping row for now, and not setting 'taxpath' means it needs to
            # be filtered out later.
            return
        else:
            raise

    if current_rank == "no rank":
        # need to figure out based on parent
        parent = taxo.parent(current_taxid)
        parent_rank = taxo.rank(parent)
        if parent_rank == "species":
            # we have a strain
            row["rank"] = "strain"
        elif parent_rank == "genus":
            # it might be a species-like rank,
            # but we should leave empty for OPAL
            row["rank"] = "no rank"
        elif parent_rank == "subspecies":
            # we have a strain
            row["rank"] = "strain"
        # NTP TODO : what do we actually want forma/varietas/etc to map to? should we keep them?
        elif parent_rank == "no rank" or parent_rank == "forma" or parent_rank == "varietas":
            # seen this with taxid 1620419, which has lineage
            # no rank|no rank|subspecies|species|genus|family|order|class|phylum|superkingdom
            parent = taxo.parent(parent)
            parent_rank = taxo.rank(parent)
            if parent_rank in ("species", "subspecies"):
                row["rank"] = "strain"
            else:
                # ntp: try going up one more
                parent = taxo.parent(parent)
                parent_rank = taxo.rank(parent)
                if parent_rank in ("species", "subspecies"):
                    row["rank"] = "strain"
                else:
                    raise Exception("TODO other ranks testing for parent")
    elif current_rank == "species":
        row["rank"] = "species"
    elif current_rank == "subspecies":
        # TODO: should probably drop it from profile if subspecies not in ranks,
        # but need to do it later? (can't remove row during apply)
        row["rank"] = "no rank"
    # NTP (same todo as above): what do we actually want forma/varietas/etc to map to? should we keep them?
    elif current_rank == "forma" or current_rank == "varietas":
        #In botanical nomenclature, a form (forma, plural formae)
        # is one of the "secondary" taxonomic ranks, below that of variety,
        # which in turn is below that of species; it is an infraspecific taxon.

        #In botanical nomenclature, variety (abbreviated var.; in Latin: varietas)
        # is a taxonomic rank below that of species and subspecies, but above that of form.
        # As such, it gets a three-part infraspecific name
        row["rank"] = "no rank"
    else:
        raise Exception("TODO other ranks testing for current")

    # uses taxonomy pkg
    lineage = taxo.lineage(current_taxid)
    valid_ranks = {}

    if row["rank"] == "strain" and "strain" in ranks:
        valid_ranks["strain"] = current_taxid

    for l in lineage:
        current_rank = taxo.rank(l)
        if current_rank in ranks:
            valid_ranks[current_rank] = l

    final_lineage = []
    final_taxpath = []
    for rank in ranks:
        if rank in valid_ranks:
            name = taxo.name(valid_ranks[rank])
            lp = sourmash.lca.lca_utils.LineagePair(rank, name)
            final_lineage.append(lp)
            #if rank == "superkingdom":
            #    row["charcoal_lineage"] = f"d__{name}"
            final_taxpath.append(valid_ranks[rank])
        else:
            final_taxpath.append("")
            lp = sourmash.lca.lca_utils.LineagePair(rank, "")
            final_lineage.append(lp)

    row["taxpath"] = "|".join(final_taxpath)
    row["lineage"] = ";".join(sourmash.lca.lca_utils.zip_lineage(final_lineage,include_strain=True,truncate_empty=True))
    return row

def main(args):
    taxdump_dir= args.taxdump_dir
    ranks=args.ranks.split("|")
    taxid_col=args.taxid_column
    # load local taxdump files names.dmp and nodes.dmp into taxonomy
    taxo = taxonomy.Taxonomy.from_ncbi(os.path.join(taxdump_dir, "nodes.dmp"), os.path.join(taxdump_dir, "names.dmp"))

    #read in csv
    queryDF= pd.read_csv(args.query_csv)
    # find lineage from taxid
    queryDF = queryDF.apply(get_row_taxpath, axis=1, args=(taxo, ranks, taxid_col))

    # write full csv
    queryDF.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("query_csv")
    p.add_argument("--output-csv")
    p.add_argument("--taxdump_dir", default="ncbi_taxonomy")
    p.add_argument("--ranks", default="superkingdom|phylum|class|order|family|genus|species|strain")
    p.add_argument("--taxid-column", default="tax_id", help="ncbi taxid column")
    p.add_argument("--identifier-column", default="accession", help="column in query_csv that can be used to identify corresponding filename")
    args = p.parse_args()
    if not args.output_csv:
        args.output_csv=(args.query_csv).rsplit(".", 1)[0] + ".lineage.csv"
    sys.exit(main(args))
