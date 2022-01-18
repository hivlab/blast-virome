import os

pepfile: "config/pep.yaml"
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"
TAXON_DB = os.getenv("TAXON_DB")


def get_contigs(wildcards):
    return pep.sample_table.loc[(pep.sample_table["sample_name"]==wildcards.sample) & (pep.sample_table["workflow"]==wildcards.workflow),"contigs"]

rule all:
    input:
        expand(["output/{sample}/{workflow}/contigs.fa", "output/{sample}/{workflow}/viruses.csv"], zip, sample=pep.sample_table["sample_name"], workflow=pep.sample_table["workflow"])

rule reformat:
    input:
        get_contigs,
    output:
        "output/{sample}/{workflow}/contigs.fa",
    conda:
        f"{WRAPPER_PREFIX}/v0.2/subset_fasta/environment.yaml"
    shell:
        "python -u scripts/fix_fasta.py --input {input[0]} --output {output[0]}"

include: "rules/blast.smk"