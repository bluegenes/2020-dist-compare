import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir

out_dir = "output.shared_kmers"
logs_dir = os.path.join(out_dir, "logs")
index_dir = "/home/ntpierce/sourmash_databases/gtdb-r95/index"

lca_indices = glob.glob(os.path.join(index_dir, "*lca.json.gz"))
lca_basenames =  [os.path.basename(idx).rsplit(".lca.json.gz", 1)[0] for idx in lca_indices]


rule all:
    input: 
        expand(os.path.join(out_dir, "{basename}.shared-kmers.csv"), basename=lca_basenames),

rule count_shared_kmers:
    input: 
        os.path.join(index_dir, "{basename}.lca.json.gz")
    output:
        os.path.join(out_dir, "{basename}.shared-kmers.csv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *50000,
        runtime=1200,
    log: os.path.join(logs_dir, "{basename}.shared-kmers.log")
    benchmark: os.path.join(logs_dir, "{basename}.shared-kmers.benchmark")
    conda: "envs/sourmash-dev.yml"
    shell:
        """
        python shared-kmers.py --db {input} --csv {output} 
        """
