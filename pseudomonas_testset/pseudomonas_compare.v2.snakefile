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


## WORKING HERE

# build comparison info
compareInfo = pd.read_csv(config["comparison_info"])
compareD=defaultdict()
for comparison rank in cInfo["lowest_common_rank"].unique():
    rank_set = cInfo[cInfo["lowest_common_rank"] == comparison_rank]
    # get anchor acc
    for anchor_acc in rank_set["anchor_acc"].unique(): # need to get value
    # get comparison acc list
        comparison_list = rank_set[rank_set["anchor_acc"] == anchor_acc]["compare_accs"].to_list()
    compareD[comparison_rank][anchor_acc] = comparison_list


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
        expand(os.path.join(compare_dir, "pseudomonas.genomic.{alphak}.anchor-compare.csv.gz"), alphak=genomic_alphaksizes),
        expand(os.path.join(compare_dir, "pseudomonas.protein.{alphak}.anchor-compare.csv.gz"), alphak=protein_alphaksizes),
        os.path.join(compare_dir, "fastani-compare", "pseudomonas.genomic.fastani.tsv"),
        os.path.join(compare_dir, "compareM", "aai/aai_summary.tsv"),

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


## todo:: make these comparison-specific (for compareM and fastANI)
localrules: write_genomic_fastani_fastalist, write_genomic_compareM_fastalist, write_protein_compareM_fastalist
rule write_genomic_fastani_fastalist:
    input:
        lambda w: ancient(expand(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"), accession= pseudomonas_accessions)) 
    output: os.path.join(compare_dir, "pseudomonas.genomic.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

localrules: build_filepaths_for_fastani
rule build_filepaths_for_fastani:
	output: os.path.join(out_dir, "fastani", "{path}", "{path}.filepaths.txt")
	run:
		with open(str(output), "w") as out:
			acc_list = path2acc[wildcards.path]
			for acc in acc_list:
				fn = os.path.join(data_dir, lineages_info.at[acc, 'filename'])
				out.write(f"{fn}\n")

rule write_genomic_compareM_fastalist:
#rule nucl_build_filepaths_for_compareM:
	output:
		os.path.join(out_dir, "compareM", "pseudomonas", "species.filepaths.txt")
	params:
		outdir = lambda w: os.path.join(out_dir, "compareM", w.path)
	group: "nuclcompareM"
	run:
		with open(str(output), "w") as out:
			acc_list = path2acc[wildcards.path]
			for acc in acc_list:
				fn = os.path.join(data_dir, lineages_info.at[acc, 'filename'])
				fn_gunzip = os.path.join(params.outdir, os.path.basename(fn).rsplit(".gz")[0])
				shell("gunzip -c {fn} > {fn_gunzip}")
				out.write(f"{fn_gunzip}\n")
  #   input:
  #      lambda w: ancient(expand(os.path.join(out_dir, "data/genomic/{accession}_genomic.fna.gz"), accession= pseudomonas_accessions)) 
  #  output: os.path.join(compare_dir, "pseudomonas.genomic.fastalist.txt")
  #  run:
  #      with open(str(output), "w") as outF:
  #          for inF in input:
  #              outF.write(str(inF) + "\n")

rule write_protein_compareM_fastalist:
    input:
        lambda w: ancient(expand(os.path.join(out_dir, "data/protein/{accession}_protein.faa.gz"), accession= pseudomonas_accessions))
    output: os.path.join(compare_dir, "pseudomonas.protein.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


### TODO: Instead of writing full siglists, let's write the comparison-level siglists
localrules: write_genomic_siglists, write_protein_siglists
rule write_genomic_siglists:
    input:
        lambda w: expand(os.path.join(out_dir, "genomic", "signatures", "{accession}.genomic.sig"), accession = pseudomonas_accessions), 
    output: os.path.join(compare_dir, "pseudomonas.genomic.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


rule write_protein_siglists:
    input:
        lambda w: expand(os.path.join(out_dir, "protein", "signatures", "{accession}.protein.sig"), accession = pseudomonas_accessions),
    output: os.path.join(compare_dir, "{dataset}.protein.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


# TODO: Make these comparison-specific!
rule taxon_compare_genomic:
    input:
        siglist = os.path.join(compare_dir, "pseudomonas.genomic.siglist.txt")
    output:
        os.path.join(compare_dir, "pseudomonas.genomic.{alphabet}-k{ksize}.anchor-compare.csv.gz"),
    params:
        taxinfo_csv = pseudomonas_csv,
        sigdir= os.path.join(out_dir, "genomic", "signatures"),
        scaled= lambda w: min(alphabet_info[w.alphabet]["scaled"]),
        sig_ext = ".genomic.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare" + expt, "pseudomonas.genomic.{alphabet}-k{ksize}.anchor-compare.log")
    benchmark: os.path.join(logs_dir, "compare" + expt, "pseudomonas.genomic.{alphabet}-k{ksize}.anchor-compare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        python compare-pseudomonas.py --taxinfo-csv {params.taxinfo_csv} \
               --compare-alphabet {wildcards.alphabet} --ksize {wildcards.ksize} \
               --scaled {params.scaled} --sigdir {params.sigdir} \
               --sig-ext {params.sig_ext} --output-csv {output} > {log} 2>&1
        """

rule taxon_compare_protein:
    input:
        siglist = os.path.join(compare_dir, "pseudomonas.protein.siglist.txt")
    output:
        os.path.join(compare_dir, "pseudomonas.protein.{alphabet}-k{ksize}.anchor-compare.csv.gz"),
    params:
        taxinfo_csv = pseudomonas_csv,
        sigdir= os.path.join(out_dir, "protein", "signatures"),
        scaled= lambda w: min(alphabet_info[w.alphabet]["scaled"]),
        sig_ext = ".protein.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare"+expt, "pseudomonas.protein.{alphabet}-k{ksize}.anchor-compare.log")
    benchmark: os.path.join(logs_dir, "compare" +expt, "pseudomonas.protein.{alphabet}-k{ksize}.anchor-compare.benchmark")
    conda: "/home/ntpierce/2020-distance-compare/envs/sourmash-dev.yml"
    shell:
        """
        python compare-pseudomonas.py  --taxinfo-csv {params.taxinfo_csv} \
               --compare-alphabet {wildcards.alphabet} --ksize {wildcards.ksize} \
               --scaled {params.scaled} --sigdir {params.sigdir} \
               --sig-ext {params.sig_ext} --output-csv {output} > {log} 2>&1
        """


# prodigal cant use gzipped nucl files!!(???)
rule nucl_AAI_via_compareM:
	input:
		os.path.join(out_dir, "compareM", "{path}/{path}.filepaths.txt")
	output:
		os.path.join(out_dir, "compareM", "{path}/aai/aai_summary.tsv"),
	params:
		proteins_cmd = "",
		file_ext = ".fna",
		outdir = lambda w: os.path.join(out_dir, "compareM", w.path),
		fna_filepath = lambda w: os.path.join(out_dir, "compareM", w.path, "*fna"),
		# sigh, this is hacky. should probably use snakemake temp files
	group: "nuclcompareM"
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *5000,
		runtime=60,
	log: os.path.join(logs_dir, "compareM/paths", "{path}.compareM.log")
	benchmark: os.path.join(logs_dir, "compareM/paths", "{path}.compareM.benchmark")
	conda: "envs/compareM-env.yml"
	shell:
		"""
		comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
		rm -rf {params.fna_filepath}
		"""


def get_genome_info(w):
	anchor_acc = pathinfo[(pathinfo["path"] == w.path) & (pathinfo["rank"] == "species")].index[0]
	anchor_g = os.path.join(data_dir, lineages_info.at[anchor_acc, 'filename'])
	path_glist =  os.path.join(out_dir, "fastani", f"{w.path}/{w.path}.filepaths.txt")
	return {"anchor_genome": anchor_g, "path_genomes": path_glist}

rule compare_via_fastANI:
	input:
		unpack(get_genome_info)
	output: os.path.join(out_dir, "fastani", "{path}/{path}.fastani.tsv"),
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *5000,
		runtime=1200,
	log: os.path.join(logs_dir, "fastani", "{path}/{path}.fastani.log")
	benchmark: os.path.join(logs_dir, "fastani", "{path}/{path}.fastani.benchmark")
	conda: "envs/fastani-env.yml"
	shell:
		"""
		fastANI -q {input.anchor_genome:q} --rl {input.path_genomes:q} -o {output} > {log} 2>&1
		"""

localrules: write_fastani_result_csv
rule write_fastani_result_csv:
	input: expand(os.path.join(out_dir, "fastani", "{path}/{path}.fastani.tsv"), path=path_names)
	output: os.path.join(out_dir, "fastani", "{basename}.path-fastani.filecsv"),
	run:
		with open(str(output), "w") as out:
			for path in path_names:
				out.write(f"{path},{out_dir}/fastani/{path}/{path}.fastani.tsv\n")

localrules: aggregate_fastani_results
rule aggregate_fastani_results:
	input:
		fastani=os.path.join(out_dir, "fastani", "{basename}.path-fastani.filecsv"),
		paths=config["evolpaths"],
	output: os.path.join(out_dir, "fastani", "{basename}.path-fastani.csv.gz"),
	log: os.path.join(logs_dir, "fastani", "{basename}.path-fastani.aggregate.log")
	benchmark: os.path.join(logs_dir, "fastani", "{basename}.path-fastani.aggregate.benchmark")
	shell:
		"""
		python aggregate-fastani-results.py --fastani-filecsv {input.fastani} \
											 --path-info {input.paths} \
											 --output-csv {output} > {log} 2>&1
		"""

# compareM proteins CAN handle gzipped files (??)
localrules: build_filepaths_for_compareM
rule protein_build_filepaths_for_compareM:
	output:
		os.path.join(out_dir, "compareM", "{path}", "{path}.filepaths.txt")
	run:
		with open(str(output), "w") as out:
			acc_list = path2acc[wildcards.path]
			for acc in acc_list:
				fn = os.path.join(data_dir, lineages_info.at[acc, 'filename'])
				out.write(f"{fn}\n")


rule protein_AAI_via_compareM:
	input:
		os.path.join(out_dir, "compareM", "{path}/{path}.filepaths.txt")
	output:
		os.path.join(out_dir, "compareM", "{path}/aai/aai_summary.tsv"),
	params:
		proteins_cmd = "--proteins" if input_type == "protein" else "",
		file_ext = ".faa.gz" if input_type == "protein" else ".fna",
		outdir = lambda w: os.path.join(out_dir, "compareM", w.path),
	threads: 1
	resources:
		mem_mb=lambda wildcards, attempt: attempt *5000,
		runtime=60,
	log: os.path.join(logs_dir, "compareM/paths", "{path}.compareM.log")
	benchmark: os.path.join(logs_dir, "compareM/paths", "{path}.compareM.benchmark")
	conda: "envs/compareM-env.yml"
	shell:
		"""
		comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
		"""


# make these work with the above
localrules: write_compareM_result_csv
rule write_compareM_result_csv:
	input:
		expand(os.path.join(out_dir, "compareM", "{path}/aai/aai_summary.tsv"), path=path_names)
	output:
		os.path.join(out_dir, "compareM", "{basename}.path-compareM.filecsv"),
	run:
		with open(str(output), "w") as out:
			for path in path_names:
				out.write(f"{path},{out_dir}/compareM/{path}/aai/aai_summary.tsv\n")

localrules: aggregate_compareM_results
rule aggregate_compareM_results:
	input:
		compareM=os.path.join(out_dir, "compareM", "{basename}.path-compareM.filecsv"),
		paths=config["evolpaths"],
	output: os.path.join(out_dir, "compareM", "{basename}.path-compareM.csv.gz"),
	log: os.path.join(logs_dir, "compareM", "{basename}.path-compareM.aggregate.log")
	benchmark: os.path.join(logs_dir, "compareM", "{basename}.path-compareM.aggregate.benchmark")
	shell:
		"""
		python aggregate-compareM-results.py --comparem-tsv-filecsv {input.compareM} \
											 --path-info {input.paths} \
											 --output-csv {output} > {log} 2>&1
		"""


