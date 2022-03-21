RANKS_OF_INTEREST = ["superkingdom", "order", "family", "genus", "species"]
VIRUSES_TAXID = 10239


# Creates to required outputs viruses.taxids and negative.taxids.
# Output directory can be changed.
# Additional negative taxids (all listed taxids except viruses) can be added via params.
# Shadow=full ensures that only required outputs will be saved.
rule get_virus_taxids:
    output:
        f"results/{VIRUSES_TAXID}.taxids",
    params:
        taxid=VIRUSES_TAXID,
    conda:
        f"{WRAPPER_PREFIX}/v0.2/blast/query/environment.yaml"
    resources:
        runtime=120,
    shell:
        "(get_species_taxids.sh -t {params.taxid} > {output} || true)"


rule megablast_virus:
    input:
        query=rules.repeatmasker.output.masked,
        taxidlist=rules.get_virus_taxids.output[0],
    output:
        out=temp("results/{workflow}/{sample}/megablast-virus.tsv"),
    params:
        program="blastn",
        task="megablast",
        db="nt_v5",
        word_size=16,
        evalue=1e-6,
        outfmt="'6 qseqid sacc staxid pident length qstart qend sstart send evalue'",
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=26000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/query"


# Filter blastn hits for the cutoff value.
rule parse_megablast_virus:
    input:
        query=rules.megablast_virus.input.query,
        blast_result=rules.megablast_virus.output.out,
    output:
        mapped=temp("results/{workflow}/{sample}/megablast-virus_mapped.tsv"),
        unmapped=temp("results/{workflow}/{sample}/megablast-virus_unmapped.fa"),
    params:
        e_cutoff=1e-6,
        outfmt=rules.megablast_virus.params.outfmt,
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/parse"


# Blastn, megablast and blastx input, output, and params keys must match commandline blast option names.
# Please see https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a for all available options.
# Blast against nt virus database.
rule blastn_virus:
    input:
        query=rules.parse_megablast_virus.output.unmapped,
        taxidlist=rules.get_virus_taxids.output[0],
    output:
        out=temp("results/{workflow}/{sample}/blastn-virus.tsv"),
    params:
        program="blastn",
        db="nt_v5",
        word_size=11,
        evalue=1e-6,
        outfmt=rules.megablast_virus.params.outfmt,
    threads: 4
    resources:
        runtime=lambda wildcards, attempt: 840 + (attempt * 120),
        mem_mb=26000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/query"


# Filter blastn hits for the cutoff value.
rule parse_blastn_virus:
    input:
        query=rules.blastn_virus.input.query,
        blast_result=rules.blastn_virus.output.out,
    output:
        mapped=temp("results/{workflow}/{sample}/blastn-virus_mapped.tsv"),
        unmapped=temp("results/{workflow}/{sample}/blastn-virus_unmapped.fa"),
    params:
        e_cutoff=1e-6,
        outfmt=rules.megablast_virus.params.outfmt,
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/parse"


rule megablast_nt:
    input:
        query=rules.parse_blastn_virus.output.unmapped,
    output:
        out=temp("results/{workflow}/{sample}/megablast-nt.tsv"),
    params:
        program="blastn",
        task="megablast",
        db="nt_v5",
        word_size=16,
        evalue=1e-6,
        outfmt=rules.megablast_virus.params.outfmt,
    threads: 8
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
        mem_mb=96000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/query"


# Filter blastn hits for the cutoff value.
rule parse_megablast_nt:
    input:
        query=rules.megablast_nt.input.query,
        blast_result=rules.megablast_nt.output.out,
    output:
        mapped=temp("results/{workflow}/{sample}/megablast-nt_mapped.tsv"),
        unmapped=temp("results/{workflow}/{sample}/megablast-nt_unmapped.fa"),
    params:
        e_cutoff=1e-6,
        outfmt=rules.megablast_virus.params.outfmt,
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/parse"


# Filter sequences by division id.
# Saves hits with division id
BLAST = ["megablast-virus", "blastn-virus", "megablast-nt"]


rule classify_all:
    input:
        expand(
            "results/{{workflow}}/{{sample}}/{blastresult}_mapped.tsv", blastresult=BLAST
        ),
    output:
        temp("results/{workflow}/{sample}/all.csv"),
    params:
        pp_sway=1,
        ranks_of_interest=RANKS_OF_INTEREST,
        dbfile=TAXON_DB,
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=8000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/blast/taxonomy"


# Split classification results into viruses and non-viral
rule filter_viruses:
    input:
        "results/{workflow}/{sample}/all.csv",
    output:
        viral="results/{workflow}/{sample}/viruses.csv",
        non_viral="results/{workflow}/{sample}/non-viral.csv",
    params:
        ranks=RANKS_OF_INTEREST,
    resources:
        runtime=lambda wildcards, attempt: attempt * 120,
        mem_mb=8000,
    run:
        tab = concatenate_tables(input, sep=",", cols_to_integer=params.ranks)
        mask = tab.superkingdom == VIRUSES_TAXID
        mask = mask.fillna(False)
        vir = tab[mask]
        non_vir = tab[~mask]
        vir.to_csv(output.viral, index=False)
        non_vir.to_csv(output.non_viral, index=False)


# Merge unassigned sequences
rule merge_unassigned:
    input:
        expand(
            "results/{{workflow}}/{{sample}}/{blastresult}_unmapped.fa",
            blastresult=BLAST
        ),
    output:
        "results/{workflow}/{sample}/unassigned.fa",
    shell:
        "cat {input} > {output}"
