# Tantan mask of low complexity DNA sequences
rule tantan:
    input:
        rules.reformat.output[0],
    output:
        temp("results/{workflow}/{sample}/repeatmasker.fa"),
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
        masked=temp("results/{workflow}/{sample}/repeatmasker.fa.masked"),
        out=temp("results/{workflow}/{sample}/repeatmasker.fa.out"),
        tbl="results/{workflow}/{sample}/repeatmasker.fa.tbl",
    shadow: "minimal"
    params:
        extra="",
        outdir=lambda wildcards, output: os.path.dirname(output.masked),
    threads: 8
    resources:
        runtime=1440,
        mem_mb=16000,
    container:
        "docker://taavipall/repeatmasker-image"
    shell:
        """
        RepeatMasker {params.extra} -pa {threads} {input} -dir {params.outdir}
        if head -n 1 {output.out} | grep -q 'There were no repetitive sequences detected'; then
            ln -sr {input} {output.masked} && touch {output.tbl}
        fi
        """

