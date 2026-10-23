# RNA 5' ends: one count per read pair at read 2's 5' end, on read 2's strand.
# Raw counts; minus-strand values stay positive. Libraries get their own track,
# experiments pool their libraries.


rule bigwig:
    input:
        bams=lambda w: [f"{out('dedup', lib)}.bam" for lib in libraries_in(w.unit)],
        sizes=f"{REFDIR}/signal.chrom.sizes",
    output:
        plus=f"{OUTDIR}/{{experiment}}/bigwig/{{unit}}.+.bw",
        minus=f"{OUTDIR}/{{experiment}}/bigwig/{{unit}}.-.bw",
    params:
        name=lambda w: out("bigwig", w.unit),
    log:
        f"{LOGDIR}/{{experiment}}/bigwig/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/{{experiment}}/bigwig/{{unit}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        OPENBLAS_NUM_THREADS=1 figwig bam2bw {input.bams:q} \
            -s {input.sizes:q} \
            -n {params.name:q} \
            -mp \
            --rna5 read2 \
            -p {threads:q}

        sync {output:q}
        """
