"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  pseudomonas-compare.snakefile -k
"""
# download all pseudomonas genomes, do containment analysis at DNA, protein levels. 
#    - calculate mean, median sub-species, species-level, genus-level containment
# find genomes matching higher taxonomic ranks, do containment analysis
#    - calculate mean, median family, class, order level containment
## this is similar, but smaller than the full-gtdb version

import os
import re
import pandas as pd

#configfile: "conf.yml"
configfile: "conf_full.yml"
basename = "pseudomonas"
out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
expt = config.get("experiment", "")
if expt:
    expt = "_" + expt
compare_dir = os.path.join(out_dir, "compare" + expt)

pseudomonas_csv = config["accession_csv"]
pseudomonas_info = pd.read_csv(pseudomonas_csv)
pseudomonas_accessions = pseudomonas_info["accession"].tolist()
pseudomonas_info.set_index("accession", inplace=True)

compareInfo = pd.read_csv(config["comparison_info"]) 
# aggregate long form dataframe to list of accessions to compare
compareInfo = compareInfo.groupby(["anchor_acc", "lowest_common_rank", "anchor_sciname"]).agg({"compare_accs":lambda x: list(x)})#.reset_index()
#compareInfo.set_index("anchor_acc", inplace=True)# set index to anchor acc 
# access info: compareInfo.loc[("GCA_000009225.1", "species")][ "compare_accs"].values[0]
# get correct alphabets, ksizes for comparisons
alphabet_info = config["alphabet_info"]
genomic_alphaksizes, protein_alphaksizes = [],[]
for alpha in alphabet_info:
    ak = expand("{alpha}-k{k}", alpha=alpha, k=alphabet_info[alpha]["ksizes"])
    genomic_alphaksizes+=ak
    if alpha != "nucleotide":
        protein_alphaksizes+=ak

rule all:
    input: 
        # fastani
        #os.path.join(compare_dir, "fastani", "pseudomonas.fastani.csv.gz"),
        # compareM
        expand(os.path.join(compare_dir, "compareM", "pseudomonas.{input_type}.compareM.csv.gz"), input_type=["genomic", "protein"]),
        # sourmash
        expand(os.path.join(compare_dir, "taxon-compare", "{basename}.{input_type}.taxoncompare.csv.gz"), basename=basename, input_type=["genomic", "protein"]),
        
        #expand(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"), accession = pseudomonas_accessions),
        #expand(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"), accession = pseudomonas_accessions),
        #expand(os.path.join(compare_dir, "pseudomonas.genomic.{alphak}.compare.csv"), alphak=genomic_alphaksizes),
        #expand(os.path.join(compare_dir, "pseudomonas.protein.{alphak}.compare.csv"), alphak=protein_alphaksizes),
        #expand(os.path.join(compare_dir, "pseudomonas.genomic.{alphak}.anchor-compare.csv.gz"), alphak=genomic_alphaksizes),
        #expand(os.path.join(compare_dir, "pseudomonas.protein.{alphak}.anchor-compare.csv.gz"), alphak=protein_alphaksizes),
        #os.path.join(compare_dir, "fastani-compare", "pseudomonas.genomic.fastani.tsv"),
        #os.path.join(compare_dir, "compareM", "aai/aai_summary.tsv"),

rule download_ncbi_datasets_tool:
    output: "scripts/ncbi-datasets"
    shell:
        """
        # linux version
        wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets
        chmod +x {output}
        """

rule ncbi_datasets_download:
    input: "scripts/ncbi-datasets"
    output:
        genomic=protected(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz")),
        protein=protected(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"))
    params:
        tmp= lambda w: os.path.join(out_dir, w.accession) + ".zip"
    threads: 1
    resources:
        mem_mb= 2000,
        runtime= 60
    shell:
        # ~ 350 or so didn't have genomic.fna but do have unplaced.scaf.fna
        # can use *fna to get this, but I'lm a bit worried about missing some _genomic - do folders ever have both?
        # hacky solution used 1/26/2021:
        #   1. ran this with *_genomic.fna pattern first, x2 to make sure all with _genomic were copied
        #   2. then I changed *_genomic.fna to *fna and ran again to get some sequence for the remaining genomes
        #   3. reverted rule to *_genomic.fna
        """
        scripts/ncbi-datasets download genome accession {wildcards.accession} --exclude-rna --exclude-gff3 -f {params.tmp}
        unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*._genomic.fna | gzip -9 > {output.genomic}
        unzip -p {params.tmp} ncbi_dataset/data/{wildcards.accession}/*.faa | gzip -9 > {output.protein}
        rm -rf {params.tmp} 
        """

def build_sketch_params(output_type):
    sketch_cmd = ""
    if output_type == "nucleotide":
        ksizes = config["alphabet_info"]["nucleotide"]["ksizes"]
        scaled = min(config["alphabet_info"]["nucleotide"]["scaled"])
        sketch_cmd = " -p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
        return sketch_cmd
    else:
        for alpha in ["protein", "dayhoff", "hp"]:
            if alpha in config["alphabet_info"].keys():
                ksizes = config["alphabet_info"][alpha]["ksizes"]
                scaled = min(config["alphabet_info"][alpha]["scaled"])
            sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

rule sourmash_sketch_genomic:
    input:
        ancient(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"))
    output:
        os.path.join(out_dir, "genomic", "signatures", "{accession}.genomic.sig"),
    params:
        nucl_sketch_params = build_sketch_params("nucleotide"),
        translate_sketch_params = build_sketch_params("protein"),
        nucl_sketch=os.path.join(out_dir, "genomic/signatures", "{accession}.nucleotide.sig"),
        prot_sketch=os.path.join(out_dir, "genomic/signatures", "{accession}.translate.sig"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{accession}_genomic.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{accession}_genomic.sketch.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch dna {params.nucl_sketch_params} --name {wildcards.accession:q} -o {params.nucl_sketch} {input}  2> {log}
        sourmash sketch translate {params.translate_sketch_params} --name {wildcards.accession:q} -o {params.prot_sketch} {input}  2> {log}
        sourmash sig cat {params.nucl_sketch} {params.prot_sketch} -o {output} 2>> {log}
        rm {params.nucl_sketch}
        rm {params.prot_sketch}
        """

rule sourmash_sketch_protein_input:
    input:
        ancient(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"))
    output:
        os.path.join(out_dir, "protein", "signatures", "{accession}.protein.sig"),
    params:
        sketch_params = build_sketch_params("protein"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{accession}_protein.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{accession}_protein.sketch.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch protein {params.sketch_params} --name {wildcards.accession:q} -o {output} {input} 2> {log}
        """

#####################
# fastANI comparisons
######################
## access 
#compareInfo.loc[("GCA_000009225.1", "species")][ "compare_accs"].values[0]

def get_fastani_comparison_genome_files(w):
    compare_accs = compareInfo.loc[(w.anchor, w.lcrank)]["compare_accs"].values[0]
    genome_paths = []
    for acc in compare_accs:
        genome_paths += [os.path.join(out_dir, f"data/genomic/{acc}_genomic.fna.gz")]
    return genome_paths 


localrules: write_genomic_fastani_fastalist
rule write_genomic_fastani_fastalist:
    input: ancient(get_fastani_comparison_genome_files)
    output: os.path.join(compare_dir, "fastani", "{lcrank}-anchor{anchor}", "{lcrank}-anchor{anchor}.genomic.fastalist")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

def get_fastani_comparison_info(w):
    compare_accs = compareInfo.loc[(w.anchor, w.lcrank)]["compare_accs"].values[0]
    genome_paths = []
    anchor_g = os.path.join(out_dir, f"data/genomic/{w.anchor}_genomic.fna.gz")
    c_filelist = os.path.join(compare_dir, "fastani", f"{w.lcrank}-anchor{w.anchor}", f"{w.lcrank}-anchor{w.anchor}.genomic.fastalist")
    return {"anchor_genome" : anchor_g, "comparison_filelist": c_filelist}

rule compare_via_fastANI:
    input:
        unpack(get_fastani_comparison_info)
    output: 
        os.path.join(compare_dir, "fastani", "{lcrank}-anchor{anchor}.fastani.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    log: os.path.join(logs_dir, "fastani", "{lcrank}-anchor{anchor}.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "{lcrank}-anchor{anchor}.fastani.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/fastani-env.yml"
    shell:
        """
        fastANI -q {input.anchor_genome:q} --rl {input.comparison_filelist:q} -o {output} > {log} 2>&1
        """

## aggreagate fastani results
def get_all_fastani(w):
    cInf = []
    for (anchor_acc, lcr, anchor_sciname) in compareInfo.index:
        cInf += [f"{lcr}-anchor{anchor_acc}"]
    fastani_files =  expand(os.path.join(compare_dir, "fastani", "{comparison_info}.fastani.tsv"), comparison_info =cInf)
    return fastani_files

localrules: write_fastani_result_csv
rule write_fastani_result_csv:
    input: get_all_fastani 
    output: os.path.join(compare_dir, "fastani", "pseudomonas.fastani.filecsv"),
    run:
        with open(str(output), "w") as out:
            for inF in input:
                comparison_name = os.path.basename(str(inF)).rsplit(".fastani.tsv")[0]
                out.write(f"{comparison_name},{str(inF)}\n")

localrules: aggregate_fastani_results
rule aggregate_fastani_results:
    input:
        fastani=os.path.join(compare_dir, "fastani", "pseudomonas.fastani.filecsv"),
        comparison_info=config["comparison_info"],
    output: os.path.join(compare_dir, "fastani", "pseudomonas.fastani.csv.gz"),
    log: os.path.join(logs_dir, "fastani", "pseudomonas.fastani.aggregate.log")
    benchmark: os.path.join(logs_dir, "fastani", "pseudomonas.fastani.aggregate.benchmark")
    shell:
        """
        python aggregate-taxon-fastani-results.py --fastani-filecsv {input.fastani} \
                                                  --comparison-info {input.comparison_info} \
                                                  --output-csv {output} > {log} 2>&1
        """


#####################
# compareM comparisons
######################

def get_compareM_protein_fastas(w):
    compare_accs = compareInfo.loc[(w.anchor, w.lcrank)]["compare_accs"].values[0]
    compare_fastas = []
    anchor_g = os.path.join(out_dir, f"data/protein/{w.anchor}_protein.faa.gz")
    for acc in compare_accs:
        compare_fastas += [os.path.join(out_dir, f"data/protein/{acc}_protein.faa.gz")]
    return compare_fastas + [anchor_g]

rule write_protein_compareM_fastalist:
    input: ancient(get_compareM_protein_fastas)
    output: os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}", "protein", "pseudomonas.{lcrank}-anchor{anchor}.fastalist")
    params:
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.lcrank}-anchor{w.anchor}", "protein")
    run:
        with open(str(output), "w") as out:
            for fn in input:
                out.write(f"{str(fn)}\n")

rule protein_AAI_via_compareM:
    input:
        os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}", "protein", "pseudomonas.{lcrank}-anchor{anchor}.fastalist")
    output:
        os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}/protein/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "--proteins",
        file_ext = ".faa.gz",
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.lcrank}-anchor{w.anchor}/protein"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=60,
    log: os.path.join(logs_dir, "compareM/protein", "{lcrank}-anchor{anchor}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/protein", "{lcrank}-anchor{anchor}.compareM.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/compareM-env.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
        """

## nucleotide compareM ##
def get_compareM_genome_fastas(w):
    compare_accs = compareInfo.loc[(w.anchor, w.lcrank)]["compare_accs"].values[0]
    compare_fastas = []
    anchor_g = os.path.join(out_dir, f"data/genomic/{w.anchor}_genomic.fna.gz") 
    for acc in compare_accs:
        compare_fastas += [os.path.join(out_dir, f"data/genomic/{acc}_genomic.fna.gz")]
    return compare_fastas + [anchor_g]

rule write_genomic_compareM_fastalist:
    input: ancient(get_compareM_genome_fastas)
    output: os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}", "genomic", "pseudomonas.{lcrank}-anchor{anchor}.fastalist")
    params:
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.lcrank}-anchor{w.anchor}", "genomic")
    group: "nuclcompareM"
    run:
        with open(str(output), "w") as out:
            for fn in input:
                fn_gunzip = os.path.join(params.outdir, os.path.basename(fn).rsplit(".gz")[0])
                shell("gunzip -c {fn} > {fn_gunzip}")
                out.write(f"{fn_gunzip}\n")

# prodigal cant use gzipped nucl files!!(???)
rule nucl_AAI_via_compareM:
    input:
        os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}", "genomic", "pseudomonas.{lcrank}-anchor{anchor}.fastalist")
    output:
        os.path.join(compare_dir, "compareM", "{lcrank}-anchor{anchor}/genomic/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "",
        file_ext = ".fna",
        outdir = lambda w: os.path.join(compare_dir, "compareM", f"{w.lcrank}-anchor{w.anchor}/genomic"),
        # sigh, this is hacky. should probably use snakemake temp files
        fna_filepath = lambda w: os.path.join(compare_dir, "compareM", f"{w.lcrank}-anchor{w.anchor}/genomic", "*fna"),
    group: "nuclcompareM"
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=60,
    log: os.path.join(logs_dir, "compareM/genomic", "{lcrank}-anchor{anchor}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/genomic", "{lcrank}-anchor{anchor}.compareM.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/compareM-env.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
        rm -rf {params.fna_filepath}
        """


## aggreagate compareM results
def get_all_compareM(w):
    cInf = []
    for (anchor_acc, lcr, lct) in compareInfo.index:
        cInf += [f"{lcr}-anchor{anchor_acc}"]
    compareM_files = expand(os.path.join(compare_dir, "compareM", "{comparison_info}/{inp}/aai/aai_summary.tsv"), comparison_info =cInf, inp=w.input_type)
    return compareM_files

rule compile_compareM_resultfiles:
    input: get_all_compareM
    output: 
        os.path.join(compare_dir, "compareM", "pseudomonas.{input_type}.compareM.filecsv",),
    run:
        with open(str(output), "w") as out:
            for inF in input:
                comparison_name = os.path.basename(str(inF).rsplit(f"/{wildcards.input_type}")[0])
                out.write(f"{comparison_name},{str(inF)}\n")

localrules: aggregate_compareM_results
rule aggregate_compareM_results:
    input:
        compareM=os.path.join(compare_dir, "compareM", "pseudomonas.{input_type}.compareM.filecsv"),
        comparison_info=config["comparison_info"],
    output: os.path.join(compare_dir, "compareM", "pseudomonas.{input_type}.compareM.csv.gz"),
    log: os.path.join(logs_dir, "compareM", "pseudomonas.{input_type}.compareM.aggregate.log")
    benchmark: os.path.join(logs_dir, "compareM", "pseudomonas.{input_type}.compareM.aggregate.benchmark")
    shell:
        """
        python aggregate-taxon-compareM-results.py --comparem-tsv-filecsv {input.compareM} \
                                                   --comparison-info {input.comparison_info} \
                                                   --output-csv {output} > {log} 2>&1
        """


#####################
# sourmash comparisons
######################

rule write_genomic_siglist:
    input:
        lambda w: expand(os.path.join(out_dir, "genomic", "signatures", "{accession}.genomic.sig"), accession = pseudomonas_accessions),
    output: os.path.join(out_dir, "pseudomonas.genomic.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule write_protein_siglist:
    input:
        lambda w: expand(os.path.join(out_dir, "protein", "signatures", "{accession}.protein.sig"), accession = pseudomonas_accessions),
    output: os.path.join(out_dir, "{dataset}.protein.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


## compare to anchor sigs ##
alpha_to_moltype = {"nucleotide": "DNA", "protein": "protein", "dayhoff": "dayhoff", "hp": "hp"}
rule taxon_compare_genomic:
    input:
        comparison_csv=config["comparison_info"],
        #lineages=config["lineages_csv"],
        sigfile=os.path.join(out_dir, "{basename}.{input_type}.siglist.txt")
    output:
        csv=os.path.join(compare_dir, "taxon-compare/{input_type}", "{basename}.{alphabet}-k{ksize}.taxoncompare.csv.gz"),
    params:
        sigdir = os.path.join(out_dir, "signatures"),
        moltype = lambda w: alpha_to_moltype[w.alphabet],
        sigext = lambda w: f".{w.input_type}.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    log: os.path.join(logs_dir, "taxon-compare/{input_type}", "{basename}.{alphabet}-k{ksize}.taxoncompare.log")
    benchmark: os.path.join(logs_dir, "taxon-compare/{input_type}", "{basename}.{alphabet}-k{ksize}.taxoncompare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/pathcompare.yml"
    shell:
        """
        python taxon-compare.py --comparison-csv {input.comparison_csv} \
        --alphabet {params.moltype} --ksize {wildcards.ksize} --sigdir {params.sigdir} \
        --siglist {input.sigfile} --sig-extension {params.sigext} --output-csv {output.csv} > {log} 2>&1
        """

localrules: aggregate_genomic_taxoncompare
rule aggregate_genomic_taxoncompare:
    input:
        expand(os.path.join(compare_dir, "taxon-compare/genomic", "{basename}.{alphak}.taxoncompare.csv.gz"), basename=basename, alphak=genomic_alphaksizes)
    output:
        os.path.join(compare_dir, "taxon-compare", "{basename}.genomic.taxoncompare.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

localrules: aggregate_protein_taxoncompare
rule aggregate_protein_taxoncompare:
    input:
        expand(os.path.join(compare_dir, "taxon-compare/protein", "{basename}.{alphak}.taxoncompare.csv.gz"), basename=basename, alphak=protein_alphaksizes)
    output:
        os.path.join(compare_dir, "taxon-compare", "{basename}.protein.taxoncompare.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

