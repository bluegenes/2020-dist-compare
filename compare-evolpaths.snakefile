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
    path2acc = paths.groupby('path')['accession'].apply(list).to_dict()
    paths.set_index("accession", inplace=True)
    return paths, path2acc



lineages_info = read_lineages(config["lineages_csv"], data_dir)
pathinfo, path2acc = read_paths(config["evolpaths"])
path_names = path2acc.keys()
sample_names = pathinfo.index.tolist()

onstart:
    print("------------------------------")
    print("Estimate similarity for 'evolutionary paths' genomes, proteomes")
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
    conditional_outputs += expand(os.path.join(out_dir, "fastani", "{basename}.path-fastani.csv.gz"),basename=basename)
    #conditional_outputs += expand(os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"), basename=basename)
 #   conditional_outputs += [os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")]
#elif input_type == "protein":
    #conditional_outputs += [os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")]


rule all:
    input: 
        conditional_outputs,
        expand(os.path.join(out_dir, "{basename}.signatures.txt"), basename=basename),
        #expand(os.path.join(out_dir, "path-compare", "{basename}.{alphak}.pathcompare.csv.gz"), basename=basename, alphak=alpha_ksizes)
        #expand(os.path.join(out_dir, "anchor-compare", "{basename}.{alphak}.anchor_containment.csv"), basename=basename, alphak=alpha_ksizes)
        expand(os.path.join(out_dir, "path-compare", "{basename}.pathcompare.csv.gz"), basename=basename),
        expand( os.path.join(out_dir, "compareM", "{basename}.path-compareM.csv.gz"), basename=basename),


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

localrules: aggregate_pathcompare
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
        

localrules: build_filepaths_for_compareM
rule build_filepaths_for_compareM:
    output:
        os.path.join(out_dir, "compareM", "{path}", "{path}.filepaths.txt")
    run:
        with open(str(output), "w") as out:
            acc_list = path2acc[wildcards.path]
            for acc in acc_list:
                fn = os.path.join(data_dir, lineages_info.at[acc, 'filename'])
                out.write(f"{fn}\n")

if input_type == "nucleotide":
    
    localrules: build_filepaths_for_fastani
    rule build_filepaths_for_fastani:
        output: os.path.join(out_dir, "fastani", "{path}", "{path}.filepaths.txt")
        run:
            with open(str(output), "w") as out:
                acc_list = path2acc[wildcards.path]
                for acc in acc_list:
                    fn = os.path.join(data_dir, lineages_info.at[acc, 'filename'])
                    out.write(f"{fn}\n")

    
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
     
rule AAI_via_compareM:
    input: 
        os.path.join(out_dir, "compareM", "{path}/{path}.filepaths.txt")
    output:
        os.path.join(out_dir, "compareM", "{path}/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "--proteins" if input_type == "protein" else "",
        file_ext = ".faa.gz" if input_type == "protein" else ".fna.gz",
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



