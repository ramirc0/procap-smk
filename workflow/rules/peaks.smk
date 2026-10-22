# Per experiment: PINTS on the pooled track, as ENCODE called peaks on its
# summed replicate bigWigs.


rule pints:
    input:
        plus=f"{OUTDIR}/{{experiment}}/{{experiment}}.+.bw",
        minus=f"{OUTDIR}/{{experiment}}/{{experiment}}.-.bw",
    output:
        multiext(
            f"{OUTDIR}/{{experiment}}/{{experiment}}_1_",
            "bidirectional_peaks.bed",
            "divergent_peaks.bed",
            "unidirectional_peaks.bed",
        ),
    params:
        outdir=lambda w: f"{OUTDIR}/{w.experiment}",
    log:
        f"{LOGDIR}/pints/{{experiment}}.txt",
    benchmark:
        f"{BENCHDIR}/pints/{{experiment}}.tsv",
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
            bed={wildcards.experiment:q}_1_${{kind}}_peaks.bed
            # PINTS writes no file for a peak type it found none of.
            if [[ -e "$tmp/$bed" ]]; then
                mv "$tmp/$bed" {params.outdir:q}/
            else
                : > {params.outdir:q}/"$bed"
            fi
        done

        sync {output:q}
        """
