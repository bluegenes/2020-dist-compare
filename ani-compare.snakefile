import os
import pandas as pd
#####
#
#Compare distance estimation on pairs of unaligned nucleotide sequences
####

# read in relevant tsv; get names
# choose the anchor sequence for comparison (choose fastani published anchor)
# all sequences --> sigs
# compute dists (each sig :: anchor) + store in tsv
# Plot true vs est pdist

configfile: "fastani-compare.yml"
out_dir = config["output_dir"]
logs_dir = out_dir + "/logs"

# get sample names; set anchors

datasets = config["datasets"]
dataset_samples = {}
for ds in datasets:
    ds_tsv = pd.read_csv(f"fastani-datasets/{ds}.ids.txt.gz", sep= "\t") #, names = ["name", "accession", "link"])
    #import pdb;pdb.set_trace()
    #sample_list = ds_tsv["name"].tolist()
    sample_list = ds_tsv.iloc[:,0].tolist() ### this doesn't work for D2!????? diff format?
    dataset_samples[ds] = sample_list

#import pdb;pdb.set_trace()

rule all: 
    input: 
        expand(os.path.join(out_dir, "compare", "{dataset}.{fasta_type}.compare.csv.gz"), dataset=datasets, fasta_type = ["dnainput", "prodigal"])

rule prodigal_translate:
    input: 
        ancient(os.path.join("fastani-datasets/{dataset}", "{name}.LargeContigs.fna"))
    output: 
        genes=os.path.join(out_dir, "{dataset}/prodigal", "{name}.genes.fasta"),
        proteins=os.path.join(out_dir, "{dataset}/prodigal", "{name}.proteins.fasta")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    shell:
        """
        prodigal -i {input} -o {output.genes} -a {output.proteins}
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
                ## if ksizes aren't given, sketch protein, dayhoff, hp at the ksizes from default config
                ksizes = config["alphabet_info"][alpha]["ksizes"]
                scaled = min(config["alphabet_info"][alpha]["scaled"])
            sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

rule sourmash_sketch_dnainput:
    input:
        ancient(os.path.join("fastani-datasets/{dataset}", "{name}.LargeContigs.fna"))
    output:
        os.path.join(out_dir, "{dataset}", "signatures", "{name}.dnainput.sig"),
    params:
        nucl_sketch_params = build_sketch_params("nucleotide"),
        translate_sketch_params = build_sketch_params("protein"),
        nucl_sketch=os.path.join(out_dir, "{dataset}/signatures", "{name}.nucleotide.sig"),
        prot_sketch=os.path.join(out_dir, "{dataset}/signatures", "{name}.translate.sig"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{dataset}__{name}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{dataset}__{name}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch dna {params.nucl_sketch_params} --name {wildcards.name:q} -o {params.nucl_sketch} {input}  2> {log}
        sourmash sketch translate {params.translate_sketch_params} --name {wildcards.name:q} -o {params.prot_sketch} {input}  2> {log}
        sourmash sig cat {params.nucl_sketch} {params.prot_sketch} -o {output} 2>> {log}
        rm {params.nucl_sketch}
        rm {params.prot_sketch}
        """

rule sourmash_sketch_prodigal_input:
    input:
        os.path.join(out_dir, "{dataset}", "prodigal", "{name}.proteins.fasta")
    output:
        os.path.join(out_dir, "{dataset}", "prodigal", "signatures", "{name}.prodigal.sig"),
    params:
        sketch_params = build_sketch_params("protein"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *2000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{dataset}__{name}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{dataset}__{name}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch protein {params.sketch_params} --name {wildcards.name:q} -o {output} {input} 2> {log}
        """

rule write_dnainput_siglist:
    input:
        #lambda w: ancient(expand("{{dataset}}/{name}.LargeContigs.fna"), name = dataset_samples[w.dataset]),
        lambda w: expand(os.path.join(out_dir, "{{dataset}}", "signatures", "{name}.dnainput.sig"), name= dataset_samples[w.dataset]), 
    output: os.path.join(out_dir, "compare", "{dataset}.dnainput.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_prodigal_siglist:
    input:
        #lambda w: ancient(expand(os.path.join(out_dir, "data/{{dataset}}/prodigal", "{name}.proteins.fasta"), name = dataset_samples[w.dataset])),
        lambda w: expand(os.path.join(out_dir, "{{dataset}}", "prodigal", "signatures", "{name}.prodigal.sig"), name = dataset_samples[w.dataset]) 
    output: os.path.join(out_dir, "compare", "{dataset}.prodigal.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule compare_dnainput:
    input: 
        siglist =  os.path.join(out_dir, "compare", "{dataset}.dnainput.siglist.txt")
    output: 
        os.path.join(out_dir, "compare", "{dataset}.dnainput.compare.csv.gz"),
    params:
        #fasta_dir = lambda w: w.dataset,
        anchor_sig = lambda w: anchor_info[w.dataset] + ".dnainput.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{dataset}.dna.log")
    benchmark: os.path.join(logs_dir, "compare", "{dataset}.dna.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python compare-sequences.py --siglist {input.siglist} \
               --anchor-sig {params.anchor_sig} \
               --output-csv {output}
        """


rule compare_prodigal:
    input:
        #dataset_info = "{dataset}.txt.gz",
        siglist = os.path.join(out_dir, "compare", "{dataset}.prodigal.siglist.txt")
    output:
        os.path.join(out_dir, "compare", "{dataset}.prodigal.compare.csv.gz"),
    params:
        #sig_dir = os.path.join(out_dir,"data/prodigal"),
        anchor_sig = lambda w: anchor_info[w.dataset] + ".prodigal.sig"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{dataset}.prodigal.log")
    benchmark: os.path.join(logs_dir, "compare", "{dataset}.prodigal.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python compare-sequences.py --siglist {input.siglist} \
               --fasta-alphabet "protein" --anchor-sig {params.anchor_sig} \
               --output-csv {output}
        """
