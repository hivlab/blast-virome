# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
        rules.reformat.output[0],
    output:
        temp("output/{sample}/{workflow}/repeatmasker.fa"),
    params:
        extra="-x N", # mask low complexity using N
    resources:
        runtime=120,
        mem_mb=8000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.2/tantan"


# Repeatmasker
# Outputs are generated from input file names by RepeatMasker
# must have file extension '.masked'
# If no repetitive sequences were detected symlink output to input file
rule repeatmasker:
    input:
        rules.tantan.output[0],
    output:
        masked=temp("output/{sample}/{workflow}/repeatmasker.fa.masked"),
        out=temp("output/{sample}/{workflow}/repeatmasker.fa.out"),
        cat=temp("output/{sample}/{workflow}/repeatmasker.fa.cat"),
        tbl="output/{sample}/{workflow}/repeatmasker.fa.tbl",
    params:
        outdir=lambda wildcards, output: os.path.dirname(output.masked),
    threads: 8
    resources:
        runtime=1440,
        mem_mb=16000,
    container:
        "docker://taavipall/repeatmasker-image"
    shell:
        "RepeatMasker {input} -dir {params.outdir}"

