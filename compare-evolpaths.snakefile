import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir

configfile: "evolpaths_config.yaml"

out_dir = config["output_dir"]
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")
data_dir = config['data_dir'].rstrip('/')
basename = config["basename"]

def sanitize_path(path):
    # expand `~`, get absolute path
    path = os.path.expanduser(path)
    path = os.path.abspath(path)
    return path

def read_lineages(samples_file, data_dir):
    samples = pd.read_csv(samples_file, dtype=str, sep=",", header=0)
    # if signame column was not given, it will be NaNs. Fill NA's with sample names.
    samples['signame'] = samples['signame'].fillna(samples['accession'])
    samples.set_index("accession", inplace=True)
    
    # Now, verify that all genome files exist
    data_dir = sanitize_path(data_dir)
    sample_list = samples["filename"].tolist()
    with open(os.path.join(out_dir, f"{basename}.filepaths.txt"), "w") as gf:
        for filename in sample_list:
            fullpath = os.path.join(data_dir, filename)
            if not os.path.exists(fullpath):
                print(f'** ERROR: genome file {filename} does not exist in {data_dir}')
            else:
                gf.write(fullpath + "\n")
    return samples

def read_paths(pathinfo_file):
    paths = pd.read_csv(pathinfo_file, dtype=str, sep="\t", header=0)
    paths.set_index("accession", inplace=True)
    return paths



lineages_info = read_lineages(config["lineages_csv"], data_dir)
pathinfo = read_paths(config["evolpaths"])
sample_names = pathinfo.index.tolist()

onstart:
    print("------------------------------")
    print("Compare jaccard and containment distances for related genomes or proteomes")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")

alphabet_info = config["alphabets"]
alpha_ksizes= []

wildcard_constraints:
    alphabet="\w+",
    ksize="\d+"


for alphabet, info in alphabet_info.items():
    ak = expand("{alpha}-k{ksize}", alpha=alphabet, ksize=info["ksizes"])
    alpha_ksizes.extend(ak)

input_type = config["input_type"]

conditional_outputs = []
if input_type == "nucleotide":
    conditional_outputs += expand(os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"), basename=basename)
    conditional_outputs += [os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")]
elif input_type == "protein":
    conditional_outputs += [os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")]


rule all:
    input: 
        conditional_outputs,
        expand(os.path.join(out_dir, "{basename}.signatures.txt"), basename=basename),
        #expand(os.path.join(out_dir, "path-compare", "{basename}.{alphak}.pathcompare.csv.gz"), basename=basename, alphak=alpha_ksizes)
        expand(os.path.join(out_dir, "path-compare", "{basename}.pathcompare.csv.gz"), basename=basename)
        #expand(os.path.join(out_dir, "anchor-compare", "{basename}.{alphak}.anchor_containment.csv"), basename=basename, alphak=alpha_ksizes)


## sketching rules ##

def build_sketch_params(output_type):
    sketch_cmd = ""
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = alphabet_info["nucleotide"]["ksizes"]
            scaled = alphabet_info["nucleotide"]["scaled"]
            # always track abund when sketching
            sketch_cmd = "dna -p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        # if input is protein, just build protein sketches
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        ## default build protein, dayhoff, hp sigs at the default ksizes from config
        ksizes = alphabet_info[alpha]["ksizes"]
        scaled = alphabet_info[alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}" + ",abund"
    return sketch_cmd

if config["input_type"] == "nucleotide":
    rule sourmash_sketch_nucleotide_input:
        input: lambda w: os.path.join(data_dir, lineages_info.at[w.sample, 'filename'])
        output:
            full_sketch=os.path.join(out_dir, "signatures", "{sample}.sig"),
        params:
            nucl_sketch_params = build_sketch_params("nucleotide"),
            translate_sketch_params = build_sketch_params("protein"),
            nucl_sketch=os.path.join(out_dir, "signatures", "{sample}.nucleotide.sig"),
            prot_sketch=os.path.join(out_dir, "signatures", "{sample}.translate.sig"),
            signame = lambda w: lineages_info.at[w.sample, 'signame'],
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *1000,
            runtime=1200,
        log: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.log")
        benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_nucl_input", "{sample}.sketch.benchmark")
        conda: "envs/sourmash-dev.yml"
        shell:
            """
            sourmash sketch {params.nucl_sketch_params} -o {params.nucl_sketch} --name {params.signame:q} {input}  2> {log}
            sourmash sketch {params.translate_sketch_params} -o {params.prot_sketch} --name {params.signame:q} {input}  2>> {log}
            sourmash sig cat {params.nucl_sketch} {params.prot_sketch} -o {output.full_sketch} 2>> {log}
            rm {params.nucl_sketch}
            rm {params.prot_sketch}
            """
else:
    rule sourmash_sketch_protein_input:
        input: lambda w: os.path.join(data_dir, lineages_info.at[w.sample, 'filename'])
        output:
            os.path.join(out_dir, "signatures", "{sample}.sig"),
        params:
            sketch_params = build_sketch_params("protein"),
            signame = lambda w: lineages_info.at[w.sample, 'signame'],
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *1000,
            runtime=1200,
        log: os.path. join(logs_dir, "sourmash_sketch_prot_input", "{sample}.sketch.log")
        benchmark: os.path.join(benchmarks_dir, "sourmash_sketch_prot_input", "{sample}.sketch.benchmark")
        conda: "envs/sourmash-dev.yml"
        shell:
            """
            sourmash sketch {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
            """

localrules: signames_to_file

rule signames_to_file:
    input:  expand(os.path.join(out_dir, "signatures", "{sample}.sig"), sample=sample_names),
    output: os.path.join(out_dir, "{basename}.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")


## compare anchor species ##


alpha_to_moltype = {"nucleotide": "DNA", "protein": "protein", "dayhoff": "dayhoff", "hp": "hp"}

rule compare_paths_to_anchor:
    input: 
        paths_csv=config["evolpaths"],
        #lineages=config["lineages_csv"],
        sigfile=os.path.join(out_dir, "{basename}.signatures.txt")
    output:
        csv=os.path.join(out_dir, "path-compare", "{basename}.{alphabet}-k{ksize}.pathcompare.csv.gz"),
    params:
        sigdir = os.path.join(out_dir, "signatures"),
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
        moltype = lambda w: alpha_to_moltype[w.alphabet],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    log: os.path.join(logs_dir, "path-compare", "{basename}.{alphabet}-k{ksize}.pathcompare.log")
    benchmark: os.path.join(logs_dir, "path-compare", "{basename}.{alphabet}-k{ksize}.pathcompare.benchmark")
    conda: "envs/pathcompare.yml"
    shell:
        """
        python path-compare.v2.py --paths-csv {input.paths_csv} \
        --alphabet {params.moltype} --ksize {wildcards.ksize} --sigdir {params.sigdir} \
        --siglist {input.sigfile} --output-csv {output.csv} > {log} 2>&1
        """

rule aggregate_pathcompare:
    input:
        expand(os.path.join(out_dir, "path-compare", "{basename}.{alphak}.pathcompare.csv.gz"), basename=basename, alphak=alpha_ksizes)
    output:
        os.path.join(out_dir, "path-compare", "{basename}.pathcompare.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)
        


if input_type == "nucleotide":
    rule compare_via_fastANI:
        input: 
            os.path.join(out_dir, "{basename}.filepaths.txt")
        output:
            os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"),
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: attempt *300000,
            runtime=1200,
        log: os.path.join(logs_dir, "fastani", "{basename}.fastani.log")
        benchmark: os.path.join(logs_dir, "fastani", "{basename}.fastani.benchmark")
        conda: "envs/fastani-env.yml"
        shell:
            """
            fastANI --ql {input} --rl {input} -o {output} > {log} 2>&1
            """

    rule AAI_via_compareM:
        input: 
            expand(os.path.join(out_dir, "{basename}.filepaths.txt"), basename=basename)
        output:
            os.path.join(out_dir, "compareM", "aai/aai_summary.tsv"),
        threads: 32
        resources:
            mem_mb=lambda wildcards, attempt: attempt *200000,
            runtime=1200,
        params:
            outdir = os.path.join(out_dir, "compareM")
        log: os.path.join(logs_dir, "compareM", "compareM.log")
        benchmark: os.path.join(logs_dir, "compareM", "compareM.benchmark")
        #shadow: "shallow"
        conda: "envs/compareM-env.yml"
        shell:
            """
            comparem aai_wf --cpus {threads} --file_ext ".fna.gz"  --sensitive {input} {params.outdir} > {log} 2>&1
            """
elif input_type == "protein":
    rule AAI_via_compareM:
        input: 
            expand(os.path.join(out_dir, "{basename}.filepaths.txt"), basename = basename)
        output:
            os.path.join(out_dir, "compareM", "aai/aai_summary.tsv"),
        threads: 32
        resources:
            mem_mb=lambda wildcards, attempt: attempt *200000,
            runtime=1200,
        params:
            outdir = os.path.join(out_dir, "compareM")
        log: os.path.join(logs_dir, "compareM", "compareM.log")
        benchmark: os.path.join(logs_dir, "compareM", "compareM.benchmark")
        #shadow: "shallow"
        conda: "envs/compareM-env.yml"
        shell:
            """
            comparem aai_wf --cpus {threads} --proteins --file_ext ".faa.gz" --sensitive {input} {params.outdir} > {log} 2>&1
            """


