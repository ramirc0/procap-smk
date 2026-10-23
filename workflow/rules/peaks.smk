# Per experiment: PINTS on the pooled track, as ENCODE called peaks on its
# summed replicate bigWigs.


rule pints:
    input:
        plus=f"{OUTDIR}/{{experiment}}/bigwig/{{experiment}}.+.bw",
        minus=f"{OUTDIR}/{{experiment}}/bigwig/{{experiment}}.-.bw",
    output:
        multiext(
            f"{OUTDIR}/{{experiment}}/pints/{{experiment}}",
            ".bidirectional.bed",
            ".divergent.bed",
            ".unidirectional.bed",
        ),
    params:
        prefix=lambda w: out("pints", w.experiment),
    log:
        f"{LOGDIR}/{{experiment}}/pints/{{experiment}}.txt",
    benchmark:
        f"{BENCHDIR}/{{experiment}}/pints/{{experiment}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # PINTS reuses the .npy caches an interrupted run leaves in --save-to.
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT

        pints_caller \
            --bw-pl {input.plus:q} \
            --bw-mn {input.minus:q} \
            --save-to "$tmp" \
            --file-prefix {wildcards.experiment:q} \
            --exp-type PROcap \
            --min-lengths-opposite-peaks 5 \
            --thread {threads:q} \
            --dont-check-updates
        cat "$tmp"/peakcalling_*.log
        for kind in bidirectional divergent unidirectional; do
            bed="$tmp"/{wildcards.experiment:q}_1_${{kind}}_peaks.bed
            # PINTS writes no file for a peak type it found none of.
            if [[ -e "$bed" ]]; then
                mv "$bed" {params.prefix:q}.$kind.bed
            else
                : > {params.prefix:q}.$kind.bed
            fi
        done

        sync {output:q}
        """


# The PRO-cap atlas training peaks: bidirectional and unidirectional calls in
# ProCapNet's 8-column format (chrom, start, end, strand, confidence, type,
# + summits, - summits), as procap-atlas src/preprocess/_merge_uni_bi_peaks.py.
rule merge_peaks:
    input:
        bi=f"{OUTDIR}/{{experiment}}/pints/{{experiment}}.bidirectional.bed",
        uni=f"{OUTDIR}/{{experiment}}/pints/{{experiment}}.unidirectional.bed",
    output:
        f"{OUTDIR}/{{experiment}}/merge_peaks/{{experiment}}.bed",
    log:
        f"{LOGDIR}/{{experiment}}/merge_peaks/{{experiment}}.txt",
    benchmark:
        f"{BENCHDIR}/{{experiment}}/merge_peaks/{{experiment}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        # Unidirectional rows come first so that the stable sort breaks ties as
        # the atlas does.
        {{
            awk -v OFS='\t' '{{
                print $1, $2, $3, $6, ".", "Unidirectional",
                    ($6 == "+" ? $(NF - 1) : "."), ($6 == "+" ? "." : $(NF - 1))
            }}' {input.uni:q}
            awk -v OFS='\t' '{{
                print $1, $2, $3, "Both", $4, "Bidirectional", $5, $6
            }}' {input.bi:q}
        }} | LC_ALL=C sort -s -k1,1 -k2,2n > {output:q}

        sync {output:q}
        """
