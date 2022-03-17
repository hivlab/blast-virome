__author__ = "Taavi Päll"
__copyright__ = "Copyright 2022, Hivlab"
__email__ = "taavi.pall@ut.ee"
__license__ = "MIT"

import os

pepfile: "config/pep.yaml"
WRAPPER_PREFIX = "https://raw.githubusercontent.com/avilab/virome-wrappers"
TAXON_DB = os.getenv("TAXON_DB")

# Get contig files names
def get_contigs(wildcards):
    return pep.sample_table.loc[
        (pep.sample_table["sample_name"] == wildcards.sample)
        & (pep.sample_table["workflow"] == wildcards.workflow),
        "contigs",
    ]

rule all:
    input:
        expand(["results/{workflow}/{sample}/contigs.fa", "results/{workflow}/{sample}/viruses.csv", "results/{workflow}/{sample}/unassigned.fa"], zip, sample=pep.sample_table["sample_name"], workflow=pep.sample_table["workflow"])

rule reformat:
    input:
        get_contigs,
    output:
        "results/{workflow}/{sample}/contigs.fa",
    conda:
        f"{WRAPPER_PREFIX}/v0.2/subset_fasta/environment.yaml"
    shell:
        "python -u scripts/fix_fasta.py --input {input[0]} --output {output[0]}"

include: "rules/common.smk"
include: "rules/mask.smk"
include: "rules/blast.smk"
