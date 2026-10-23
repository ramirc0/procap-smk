# Per library: merge its runs, keep unique proper pairs without a spliced mate,
# then remove UMI duplicates.


rule filter_unique:
    input:
        lambda w: [f"{out('star_align', r)}.bam" for r in RUNS_OF[w.library]],
    output:
        bam=temp(f"{OUTDIR}/{{experiment}}/filter_unique/{{library}}.bam"),
        bai=temp(f"{OUTDIR}/{{experiment}}/filter_unique/{{library}}.bam.bai"),
    log:
        f"{LOGDIR}/{{experiment}}/filter_unique/{{library}}.txt",
    benchmark:
        f"{BENCHDIR}/{{experiment}}/filter_unique/{{library}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT

        samtools merge -@ {threads:q} -u -o - {input:q} \
            | samtools view -@ {threads:q} -b -q 255 -f 3 -o "$tmp/unique.bam" -
        # A pair goes when either mate is spliced (N in its CIGAR).
        samtools view "$tmp/unique.bam" \
            | awk '$6 ~ /N/ {{ print $1 }}' \
            | sort -u \
            > "$tmp/spliced.txt"
        samtools view -@ {threads:q} -b -N "^$tmp/spliced.txt" -o {output.bam:q} "$tmp/unique.bam"
        samtools index -@ {threads:q} {output.bam:q}

        sync {output:q}
        """


rule dedup:
    input:
        bam=f"{OUTDIR}/{{experiment}}/filter_unique/{{library}}.bam",
        bai=f"{OUTDIR}/{{experiment}}/filter_unique/{{library}}.bam.bai",
    output:
        bam=f"{OUTDIR}/{{experiment}}/dedup/{{library}}.bam",
        bai=f"{OUTDIR}/{{experiment}}/dedup/{{library}}.bam.bai",
        log=f"{OUTDIR}/{{experiment}}/dedup/{{library}}.umi_tools.log",
    log:
        f"{LOGDIR}/{{experiment}}/dedup/{{library}}.txt",
    benchmark:
        f"{BENCHDIR}/{{experiment}}/dedup/{{library}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        umi_tools dedup \
            --paired \
            --umi-separator=: \
            --unpaired-reads=discard \
            --random-seed=1 \
            --log={output.log:q} \
            --stdin={input.bam:q} \
            --stdout={output.bam:q}
        samtools index {output.bam:q}

        sync {output:q}
        """
