"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  pseudomonas-containment.snakefile -k
"""
# download all pseudomonas genomes, do containment analysis at DNA, protein levels. 
#    - calculate mean, median sub-species, species-level, genus-level containment
# find genomes matching higher taxonomic ranks, do containment analysis
#    - calculate mean, median family, class, order level containment
## this is similar, but smaller than the full-gtdb version

import os
import re
import pandas as pd

configfile: "conf.yml"
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
        expand(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"), accession = pseudomonas_accessions),
        expand(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"), accession = pseudomonas_accessions),
        expand(os.path.join(compare_dir, "pseudomonas.genomic.{alphak}.compare.csv"), alphak=genomic_alphaksizes),
        expand(os.path.join(compare_dir, "pseudomonas.protein.{alphak}.compare.csv"), alphak=protein_alphaksizes),
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

localrules: write_dnainput_siglist, write_dnainput_fastalist, write_protein_siglist, write_protein_fastalist

rule write_dnainput_siglist:
    input:
        lambda w: expand(os.path.join(out_dir, "genomic", "signatures", "{accession}.genomic.sig"), accession = pseudomonas_accessions), 
    output: os.path.join(compare_dir, "pseudomonas.genomic.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_dnainput_fastalist:
    input:
        lambda w: ancient(expand(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"), accession= pseudomonas_accessions)) 
    output: os.path.join(compare_dir, "psuedomonas.genomic.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule write_protein_siglist:
    input:
        lambda w: expand(os.path.join(out_dir, "protein", "signatures", "{accession}.protein.sig"), accession = pseudomonas_accessions),
    output: os.path.join(compare_dir, "{dataset}.protein.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule write_protein_fastalist:
    input:
        lambda w: ancient(expand(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"), accession= pseudomonas_accessions))
    output: os.path.join(compare_dir, "psuedomonas.protein.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule compare_genomic:
    input: 
        siglist = os.path.join(compare_dir, "pseudomonas.genomic.siglist.txt")
    output: 
        np=os.path.join(compare_dir, "pseudomonas.genomic.{alphabet}-k{ksize}.compare.np"),
        csv=os.path.join(compare_dir, "pseudomonas.genomic.{alphabet}-k{ksize}.compare.csv"),
    params:
        alpha_cmd = lambda w: config["alphabet_info"][w.alphabet]["alpha_cmd"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare"+ expt, "pseudomonas.genomic.{alphabet}-k{ksize}.compare.log")
    benchmark: os.path.join(logs_dir, "compare"+ expt, "pseudomonas.genomic.{alphabet}-k{ksize}.compare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare --containment --from-file {input.siglist} \
                 --ksize {wildcards.ksize} {params.alpha_cmd} \
                 -o {output.np} --csv {output.csv}
        """

rule compare_protein:
    input:
        siglist = os.path.join(compare_dir, "pseudomonas.protein.siglist.txt")
    output:
        np=os.path.join(compare_dir, "pseudomonas.protein.{alphabet}-k{ksize}.compare.np"),
        csv=os.path.join(compare_dir, "pseudomonas.protein.{alphabet}-k{ksize}.compare.csv"),
    params:
        alpha_cmd = lambda w: config["alphabet_info"][w.alphabet]["alpha_cmd"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare"+ expt, "pseudomonas.protein.{alphabet}-k{ksize}.protein.log")
    benchmark: os.path.join(logs_dir, "compare"+ expt, "pseudomonas.protein.{alphabet}-k{ksize}.protein.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare --containment --from-file {input.siglist} \
                 --ksize {wildcards.ksize} {params.alpha_cmd} \
                 -o {output.np} --csv {output.csv}
        """

rule compare_via_fastANI:
    input:
        fastalist = os.path.join(compare_dir, "pseudomonas.genomic.fastalist.txt")
    output:
        os.path.join(compare_dir, "fastani-compare", "pseudomonas.genomic.fastani.tsv"),
    threads: 10
    resources:
        mem_mb=lambda wildcards, attempt: attempt *200000,
        runtime=6000,
    log: os.path.join(logs_dir, "fastani"+ expt, "pseudomonas.genomic.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani"+ expt, "pseudomonas.genomic.fastani.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/fastani-env.yml"
    shell:
        """
        fastANI -ql {input.fastalist} --rl {input.fastalist} -t {threads} -o {output} >> {log} 2>&1
        """

rule AAI_via_compareM:
    input:
        os.path.join(compare_dir, "compare", "pseudomonas.protein.fastalist.txt")
    output:
        os.path.join(compare_dir, "compareM", "aai/aai_summary.tsv"),
    threads: 20
    resources:
        mem_mb=lambda wildcards, attempt: attempt *100000,
        runtime=6000,
    params:
        outdir = os.path.join(compare_dir, "compareM")
    log: os.path.join(logs_dir, "compareM"+expt, "pseudomonas.protein.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM"+expt, "pseudomonas.protein.compareM.benchmark")
    shadow: "shallow"
    conda: "/home/ntpierce/2020-distance-compare/envs/compareM-env.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} --file_ext ".proteins.fasta" --proteins --sensitive {input} {params.outdir} >> {log} 2>&1
        """
