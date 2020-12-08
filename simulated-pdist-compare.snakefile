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
nt_frequencies = ["1", "f2", "f3"]
# evolutionary model
# model parameters (i.e. GTR: six relative rates of nucleotide substitution; GTR+Γ: six rates and one Γ shape parameter)
# nogam = GTR; gamma = GTR + Γ  
evolmodels = ["nogam", "gamma"]

simulation_info = expand("data-d{d}-{freq}-{model}", d = sim_distances, freq =nt_frequencies, model=evolmodels)

f1_info = expand("{siminfo}-seed{seed}-seq{seq}", siminfo = simulation_info, seed = f1_seqnums, seq=[1,2])
f2_info = expand("{siminfo}-seed{seed}-seq{seq}", siminfo = simulation_info, seed = f2_seqnums, seq=[1,2])
f3_info = expand("{siminfo}-seed{seed}-seq{seq}", siminfo = simulation_info, seed = f3_seqnums, seq=[1,2])
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
        expand(out_dir + "/dna-input/sigs/{siminfo}.sig", siminfo= simulated_readinfo),
        expand(out_dir + "/prodigal-input/sigs/{siminfo}.sig", siminfo= simulated_readinfo),
        #out_dir + "/results/estimated-distances.dnainput.csv",
        #out_dir + "/results/estimated-distances.prodigal.csv" ,

rule download_tsv:
    output: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    params:
        url=lambda w: f"https://zenodo.org/record/4034462/files/{w.siminfo}.tsv.xz?download=1"
    shell:
        """
        wget -O {output} {params.url}
        """

rule simreads_to_fasta_f1:
    input: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f1_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

rule simreads_to_fasta_f2:
    input: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f2_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """
rule simreads_to_fasta_f3:
    input: os.path.join(out_dir,"data/{siminfo}.tsv.xz")
    output: expand(os.path.join(out_dir, "data/simreads", "{{siminfo}}-seed{seed}-seq{seq}.fasta"), seed = f3_seqnums, seq=[1,2])
    params:
        outdir = os.path.join(out_dir, "data", "simreads")
    conda: "envs/pdist-env.yml"
    shell:
        """
        python simreads-to-fasta.py {input} --outdir {params.outdir}
        """

def build_sketch_params(output_type, input_type=config.get("input_type")):
    sketch_cmd = ""
    # if input is dna, build dna, translate sketches
    if input_type == "nucleotide":
        if output_type == "nucleotide":
            ksizes = alphabet_info["nucleotide"]["ksizes"]
            scaled = alphabet_info["nucleotide"]["scaled"]
            # don't need abund for these comparisons 
            sketch_cmd = "dna -p " + "k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}"# + ",abund"
            return sketch_cmd
        else:
            sketch_cmd = "translate "
    else:
        # if input is protein, just build protein sketches
        sketch_cmd = "protein "
    for alpha in ["protein", "dayhoff", "hp"]:
        ksizes = alphabet_info[alpha]["ksizes"]
        scaled = alphabet_info[alpha]["scaled"]
        sketch_cmd += " -p " + alpha + ",k=" + ",k=".join(map(str, ksizes)) + f",scaled={str(scaled)}"# + ",abund"
    return sketch_cmd

rule sourmash_sketch_nucleotide_input:
    input: os.path.join(out_dir, "data/simreads", "{siminfo}-seed{seed}-seq{seq}.fasta") 
    output:
        full_sketch=os.path.join(out_dir, "dna-input/sigs/{siminfo}-seed{seed}-seq{seq}.sig")
    params:
        nucl_sketch_params = build_sketch_params("nucleotide", input_type="nucleotide"),
        translate_sketch_params = build_sketch_params("protein", input_type="nucleotide"),
        nucl_sketch=os.path.join(out_dir, "sigs", "{siminfo}-seed{seed}-seq{seq}.nucleotide.sig"),
        prot_sketch=os.path.join(out_dir, "sigs", "{siminfo}-seed{seed}-seq{seq}.translate.sig"),
        signame = lambda w: f"data-{w.d}-{w.freq}-{w.model}-seed{w.seed}-seq{w.seq}",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{siminfo}-seed{seed}-seq{seq}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_nucl_input", "{siminfo}-seed{seed}-seq{seq}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch {params.nucl_sketch_params} -o {params.nucl_sketch} --name {params.signame} {input}  2> {log}
        sourmash sketch {params.translate_sketch_params} -o {params.prot_sketch} --name {params.signame} {input}  2>> {log}
        sourmash sig cat {params.nucl_sketch} {params.prot_sketch} -o {output.full_sketch} 2>> {log}
        rm {params.nucl_sketch}
        rm {params.prot_sketch}
        """

rule prodigal_translate:
    input: os.path.join(out_dir, "data/simreads", "{siminfo}-seed{seed}-seq{seq}.fasta")
    output: 
        genes=os.path.join(out_dir, "data/prodigal", "{siminfo}-seed{seed}-seq{seq}.genes.fasta"),
        proteins=os.path.join(out_dir, "data/prodigal", "{siminfo}-seed{seed}-seq{seq}.proteins.fasta")
    shell:
        """
        prodigal -i {input} -o {output.genes} -a {output.proteins}
        """

rule sourmash_sketch_protein_input:
    input: os.path.join(out_dir, "data/prodigal", "{siminfo}-seed{seed}-seq{seq}.proteins.fasta") 
    output: os.path.join(out_dir, "prodigal-input/sigs/{siminfo}-seed{seed}-seq{seq}.sig")
    params:
        sketch_params = build_sketch_params("protein", input_type="protein"),
        signame = lambda w: f"data-{w.d}-{w.freq}-{w.model}-seed{w.seed}-seq{w.seq}-prodigal",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path. join(logs_dir, "sourmash_sketch_prot_input", "{siminfo}-seed{seed}-seq{seq}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_prot_input", "{siminfo}-seed{seed}-seq{seq}.sketch.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash sketch {params.sketch_params} -o {output} --name {params.signame} {input} 2> {log}
        """

rule dnainput_signames_to_file_f1:
    input: expand(os.path.join(out_dir, "dna-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f1_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule dnainput_signames_to_file_f2:
    input: expand(os.path.join(out_dir, "dna-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f2_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule dnainput_signames_to_file_f3:
    input: expand(os.path.join(out_dir, "dna-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f3_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.dnainput.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule prodigal_signames_to_file_f1:
    input: expand(os.path.join(out_dir, "prodigal-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f1_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.prodigal.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule prodigal_signames_to_file_f2:
    input: expand(os.path.join(out_dir, "prodigal-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f2_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.prodigal.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule prodigal_signames_to_file_f3:
    input: expand(os.path.join(out_dir, "prodigal-input/sigs", "{{siminfo}}-seed{seed}-seq{seq}.sig"), seed = f3_seqnums, seq=[1,2])
    output: os.path.join(out_dir, "compare", "{siminfo}.prodigal.signatures.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule jaccard_compare_sigs_dnainput:
    input: os.path.join(out_dir, "compare", "{siminfo}.dnainput.signatures.txt")
    output: 
        csv=os.path.join(out_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.jaccard.csv"),
        np=os.path.join(out_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.jaccard.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.jaccard.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.jaccard.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule jaccard_compare_sigs_prodigal:
    input: os.path.join(out_dir, "compare", "{siminfo}.prodigal.signatures.txt")
    output:
        csv=os.path.join(out_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.csv"),
        np=os.path.join(out_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """


rule containment_compare_sigs_dnainput:
    input: os.path.join(out_dir, "compare", "{siminfo}.dnainput.signatures.txt")
    output:
        csv=os.path.join(out_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.containment.csv"),
        np=os.path.join(out_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.containment.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.containment.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.dnainput.{alphabet}-k{ksize}-scaled{scaled}.containment.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --containment --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

rule containment_compare_sigs_prodigal:
    input: os.path.join(out_dir, "compare", "{siminfo}.prodigal.signatures.txt")
    output:
        csv=os.path.join(out_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.containment.csv"),
        np=os.path.join(out_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.containment.np"),
    threads: 1
    params:
        alpha_cmd = lambda w: alphabet_info[w.alphabet]["alpha_cmd"],
        ksize = lambda w: int(w.ksize)*int(alphabet_info[w.alphabet]["ksize_multiplier"]),
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=6000,
    log: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.log")
    benchmark: os.path.join(logs_dir, "compare", "{siminfo}.prodigal.{alphabet}-k{ksize}-scaled{scaled}.jaccard.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        sourmash compare -o {output.np} --csv {output.csv} \
        --containment --ignore-abundance --ksize {params.ksize} \
        {params.alpha_cmd}  \
        --from-file {input}  2> {log}
        """

# write lists of files to aggregate
rule write_simulation_filelist:
    input: expand(out_dir + "data/{siminfo}.tsv.xz", siminfo=simulation_info)
    output: out_dir + "data/simulation-info.filelist.txt"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_jaccard_dnainput_filelist:
    input: expand(os.path.join(out_dir, "compare", "{siminfo}.dnainput.{aks}.jaccard.csv"), aks=alpha_ksize_scaled, siminfo=simulation_info)
    output: out_dir + "compare/dnainput-jaccard.filelist.txt"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_jaccard_prodigal_filelist:
    input: expand(os.path.join(out_dir, "compare", "{siminfo}.prodigal.{aks}.jaccard.csv"), aks=alpha_ksize_scaled, siminfo=simulation_info)
    output: out_dir + "compare/prodigal-jaccard.filelist.txt"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_containment_dnainput_filelist:
    input: expand(os.path.join(out_dir, "compare", "{siminfo}.dnainput.{aks}.containment.csv"), aks=alpha_ksize_scaled, siminfo=simulation_info)
    output: out_dir + "compare/dnainput-containment.filelist.txt"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

rule write_containment_prodigal_filelist:
    input: expand(os.path.join(out_dir, "compare", "{siminfo}.prodigal.{aks}.containment.csv"), aks= alpha_ksize_scaled, siminfo=simulation_info)
    output: out_dir + "compare/prodigal-containment.filelist.txt"
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                outF.write(str(inF) + "\n")

# aggregate info from multiple tsv files --> single csv file
rule aggregate_simulation_info:
    input: out_dir + "data/simulation-info.filelist.txt" 
    output: out_dir + "data/simulation-info.tsv"
    conda: "envs/pdist-env.yml"
    shell:
        """
        python aggregate-simulation-info.py {input} {output}
        """

# aggregate simulation info + sourmash compare info
rule aggregate_dnainput_compare:
    input:
        jaccard = out_dir + "compare/dnainput-jaccard.filelist.txt",
        containment = out_dir + "compare/dnainput-containment.filelist.txt",
        simulation_info =  out_dir + "data/simulation-info.tsv"
    output: out_dir + "results/estimated-distances.dnainput.csv"
    conda: "envs/pdist-env.yml"
    shell:
        """
        python aggregate-compare-info.py --jaccard-files {input.jaccard} --containment-files {input.containment} --siminfo {input.simulation_info} --output {output}
        """

rule aggregate_prodigal_compare:
    input:
        jaccard = out_dir + "compare/prodigal-jaccard.filelist.txt",
        containment = out_dir + "compare/prodigal-containment.filelist.txt",
        simulation_info =  out_dir + "data/simulation-info.tsv",
    output: out_dir + "results/estimated-distances.prodigal.csv"
    conda: "envs/pdist-env.yml"
    shell:
        """
        python aggregate-compare-info.py --jaccard-files {input.jaccard} --containment-files {input.containment} --siminfo {input.simulation_info} --output {output}
        """
