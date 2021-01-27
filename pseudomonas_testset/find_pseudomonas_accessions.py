"""
First, find genome assemblies with corresponding protein sequences using ncbi.datasets api
See example at https://github.com/ncbi/datasets/blob/master/examples/jupyter/ncbi-datasets-pylib/ncbi-datasets-assembly.ipynb
"""

import os
import re
import pandas as pd
import ncbi.datasets

AssemblyInfo = namedtuple('AssemblyInfo',
                          'accession, display_name, tax_id, parent_tax_id, sci_name, strain, title, submission_date, seq_length, assembly_level')

# ncbi datasets api instance
api_instance = ncbi.datasets.GenomeApi(ncbi.datasets.ApiClient())
tax_name = 'pseudomonas'
genome_summary = api_instance.assembly_descriptors_by_taxon(taxon=tax_name,limit='all')
print(f"Number of assemblies in the group '{tax_name}': {genome_summary.total_count}")

prot_info = []
acc_set=set()

for assembly in map(lambda d: d.assembly, genome_summary.assemblies):
    if not assembly.annotation_metadata:
        continue
    unique_acc = assembly.assembly_accession[3:]
    # want to ignore GCA, GCF. Refseq and genbank should have same assembly, right?!
    if unique_acc in acc_set:
        continue
    #print(assembly.assembly_accession)
    for fileinfo in assembly.annotation_metadata.file:
        if fileinfo.type == "PROT_FASTA":
     #       print("proteins found!")
            acc_set.add(unique_acc)
            prot_info.append(AssemblyInfo(accession=assembly.assembly_accession, \
                                          display_name=assembly.display_name, \
                                          tax_id=assembly.org.tax_id, \
                                          parent_tax_id=assembly.org.parent_tax_id, \
                                          sci_name = assembly.org.sci_name, \
                                          strain = assembly.org.strain, \
                                          title = assembly.org.title, \
                                          submission_date = assembly.submission_date, \
                                          seq_length = assembly.seq_length, \
                                          assembly_level = assembly.assembly_level))
            break

# convert namedtuple to pandas dataframe
infoDF = pd.DataFrame.from_records(prot_info, columns = AssemblyInfo._fields)

pseudomonas_accessions = infoDF.accession.tolist()
print(f"Number of unique assemblies with protein sequences in the group '{tax_name}': {len(pseudomonas_accessions)}")

infoDF.to_csv("pseudomonas_genomes.info.csv", index=False)
