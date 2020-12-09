#####
#
#Compare distance estimation on pairs of unaligned nucleotide sequences
####


# read in relevant tsv
# for each pair of sequences, estimate pdist
# keep track of seqID_params: true_pdist est_pdist
# Plot true vs est pdist

configfile: "simulated-pdist-compare.yml"
out_dir = config["output_dir"]
logs_dir = out_dir + "/logs"

def range_with_floats_list(start, stop, step):
    rangelist = []
    while stop > start:
        val = round(start, 2)
        rangelist.append(format(val, '.2f'))
        start += step
    return rangelist
 
# make variables that correspond to the simulated distances
sim_distances = range_with_floats_list(0.05, 1.00, 0.05)
f1_seqnums=[*range(1, 201)]
f2_seqnums=[*range(201, 401)]
f3_seqnums=[*range(401, 601)]

# frequencies
#f1: equal frequencies, i.e. freq(A) = freq(C) = freq(G) = freq(T) = 0.25,
#f2: GC-rich, i.e. freq(A) = 0.1, freq(C) = 0.3, freq(G) = 0.4, freq(T) = 0.2,
#f3: AT-rich, i.e. freq(A) = freq(T) = 0.4, freq(C) = freq(G) = 0.1.
nt_frequencies = ["f1", "f2", "f3"]
# evolutionary model
# model parameters (i.e. GTR: six relative rates of nucleotide substitution; GTR+Γ: six rates and one Γ shape parameter)
# nogam = GTR; gamma = GTR + Γ  
evolmodels = ["nogam", "gamma"]

simulation_info = expand("data-d{d}-{freq}-{model}", d = sim_distances, freq =nt_frequencies, model=evolmodels)

f1_info = expand("data-d{d}-f1-{model}-seed{seed}-seq{seq}", d = sim_distances, model=evolmodels, seed = f1_seqnums, seq=[1,2])
f2_info = expand("data-d{d}-f2-{model}-seed{seed}-seq{seq}", d = sim_distances, model=evolmodels, seed = f2_seqnums, seq=[1,2])
f3_info = expand("data-d{d}-f3-{model}-seed{seed}-seq{seq}", d = sim_distances, model=evolmodels, seed = f3_seqnums, seq=[1,2])
simulated_readinfo = f1_info + f2_info + f3_info

alphabet_info = config["alphabets"]
alpha_ksizes, alpha_ksize_scaled= [], []

wildcard_constraints:
    alphabet="\w+",
    ksize="\d+"


for alphabet, info in alphabet_info.items():
    ak = expand("{alpha}-k{ksize}", alpha=alphabet, ksize=info["ksizes"])
    aks = expand("{alpha}-k{ksize}-scaled{scaled}", alpha=alphabet, ksize=info["ksizes"], scaled = info["scaled"])
    alpha_ksizes.extend(ak)
    alpha_ksize_scaled.extend(aks)




rule all: 
    input: 
        expand(os.path.join(out_dir, "compare", "{siminfo}.{fasta_type}.compare.csv.gz"), siminfo=simulation_info, fasta_type = ["dnainput", "prodigal"])

rule download_tsv:
    output: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    params:
        url=lambda w: f"https://zenodo.org/record/4034462/files/{w.siminfo}.tsv.xz?download=1"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        wget -O {output} {params.url}
        """

rule simreads_to_fasta_f1:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f1_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

rule simreads_to_fasta_f2:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f2_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """
rule simreads_to_fasta_f3:
    input: ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz"))
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f3_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

rule prodigal_translate:
    input: 
        ancient(os.path.join(out_dir, "data/simreads", "{siminfo}-seed{seed}-seq{seq}.fasta"))
    output: 
        genes=os.path.join(out_dir, "data/prodigal", "{siminfo}-seed{seed}-seq{seq}.genes.fasta"),
        proteins=os.path.join(out_dir, "data/prodigal", "{siminfo}-seed{seed}-seq{seq}.proteins.fasta")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=1200,
    shell:
        """
        prodigal -i {input} -o {output.genes} -a {output.proteins}
        """

rule write_jaccard_dnainput_filelist:
    input: 
        f1 = ancient(expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f1_seqnums, seq=[1,2])),
        f2 = ancient(expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f2_seqnums, seq=[1,2])),
        f3 = ancient(expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f3_seqnums, seq=[1,2])),
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            all_inputs = input.f1 + input.f2 + input.f3
            for inF in all_inputs:
                outF.write(str(inF) + "\n")

rule write_jaccard_prodigal_filelist:
    input:
        f1 = ancient(expand(os.path.join(out_dir, "data/prodigal", "{{siminfo}}-seed{seed}-seq{seq}.proteins.fasta"), seed = f1_seqnums, seq=[1,2])),
        f2 = ancient(expand(os.path.join(out_dir, "data/prodigal", "{{siminfo}}-seed{seed}-seq{seq}.proteins.fasta"), seed = f2_seqnums, seq=[1,2])),
        f3 = ancient(expand(os.path.join(out_dir, "data/prodigal", "{{siminfo}}-seed{seed}-seq{seq}.proteins.fasta"), seed = f3_seqnums, seq=[1,2])),
    output: os.path.join(out_dir, "compare", "{siminfo}.prodigal.fastalist.txt")
    run:
        with open(str(output), "w") as outF:
            all_inputs = input.f1 + input.f2 + input.f3
            for inF in all_inputs:
                outF.write(str(inF) + "\n")

rule compare_dnainput:
    input: 
        simulation_info = ancient(os.path.join(out_dir,"data/{siminfo}.tsv.xz")),
        fasta_filelist = os.path.join(out_dir, "compare", "{siminfo}.dnainput.fastalist.txt")
    output: 
        os.path.join(out_dir, "compare", "{siminfo}.dnainput.compare.csv.gz"),
    params:
        fasta_dir = os.path.join(out_dir,"data/simreads")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python compare-paired-sequences.py --simulation-csv {input.simulation_info} \
               --fasta-dir {params.fasta_dir} --output-csv {output}
        """


rule compare_prodigal:
    input:
        simulation_info = os.path.join(out_dir,"data/{siminfo}.tsv.xz"),
        fasta_filelist = os.path.join(out_dir, "compare", "{siminfo}.prodigal.fastalist.txt")
    output:
        os.path.join(out_dir, "compare", "{siminfo}.prodigal.compare.csv.gz"),
    params:
        fasta_dir = os.path.join(out_dir,"data/prodigal")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python compare-paired-sequences.py --simulation-csv {input.simulation_info} \
               --fasta-dir {params.fasta_dir} --output-csv {output}
        """

# aggregate info from multiple tsv files --> single csv file
#rule aggregate_:
#    input: os.path.join(out_dir, "data/simulation-info.filelist.txt")
#    output: os.path.join(out_dir, "data/simulation-info.csv.gz")
#    conda: "envs/pdist-env.yml"
#    threads: 1
#    resources:
#        mem_mb=lambda wildcards, attempt: attempt *10000,
#        runtime=600,
#    shell:
#        """
#        python aggregate-simulation-info.py {input} {output}
#        """

