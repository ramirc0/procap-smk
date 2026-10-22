# Read, alignment and signal QC, joined into run-level tables and a MultiQC report.


rule fastqc:
    input:
        r1=lambda w: R1_OF[w.run] if w.stage == "raw" else f"{prefix(w.run)}_R1.fastq.gz",
        r2=lambda w: R2_OF[w.run] if w.stage == "raw" else f"{prefix(w.run)}_R2.fastq.gz",
    output:
        multiext(
            f"{OUTDIR}/{{run}}/{{run}}.{{stage}}",
            "_R1_fastqc.zip",
            "_R1_fastqc.html",
            "_R2_fastqc.zip",
            "_R2_fastqc.html",
        ),
    params:
        prefix=lambda w: f"{prefix(w.run)}.{w.stage}",
    wildcard_constraints:
        stage="raw|trimmed",
    log:
        f"{LOGDIR}/fastqc/{{run}}.{{stage}}.txt",
    benchmark:
        f"{BENCHDIR}/fastqc/{{run}}.{{stage}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT

        fastqc --quiet --threads {threads:q} --outdir "$tmp" {input.r1:q} {input.r2:q}
        # FastQC names its reports after the input file, minus .gz and .fastq/.fq.
        stem() {{ basename "$1" | sed -E 's/\.gz$//; s/\.(fastq|fq)$//'; }}
        for ext in zip html; do
            mv "$tmp/$(stem {input.r1:q})_fastqc.$ext" {params.prefix:q}_R1_fastqc.$ext
            mv "$tmp/$(stem {input.r2:q})_fastqc.$ext" {params.prefix:q}_R2_fastqc.$ext
        done

        sync {output:q}
        """


# A run's unfiltered BAM, or a library's deduplicated one.
rule samtools_stats:
    input:
        lambda w: f"{prefix(w.name)}.unfiltered.bam" if w.name in LIBRARY_OF else f"{prefix(w.name)}.bam",
    output:
        f"{OUTDIR}/{{name}}/{{name}}.samtools_stats.txt",
    wildcard_constraints:
        name=_alternatives(RUNS + LIBRARIES),
    log:
        f"{LOGDIR}/samtools_stats/{{name}}.txt",
    benchmark:
        f"{BENCHDIR}/samtools_stats/{{name}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools stats -@ {threads:q} {input:q} > {output:q}

        sync {output:q}
        """


rule idxstats:
    input:
        bam=f"{OUTDIR}/{{library}}/{{library}}.bam",
        bai=f"{OUTDIR}/{{library}}/{{library}}.bam.bai",
    output:
        f"{OUTDIR}/{{library}}/{{library}}.idxstats.tsv",
    log:
        f"{LOGDIR}/idxstats/{{library}}.txt",
    benchmark:
        f"{BENCHDIR}/idxstats/{{library}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        samtools idxstats {input.bam:q} > {output:q}

        sync {output:q}
        """


rule signal_qc:
    input:
        plus=f"{OUTDIR}/{{unit}}/{{unit}}.+.bw",
        minus=f"{OUTDIR}/{{unit}}/{{unit}}.-.bw",
        gtf=REFERENCES["gtf"],
        sizes=f"{REFDIR}/signal.chrom.sizes",
        peaks=lambda w: [] if w.unit in RUNS_OF else peak_beds(w.unit),
    output:
        f"{OUTDIR}/{{unit}}/{{unit}}.signal_qc.tsv",
    log:
        f"{LOGDIR}/signal_qc/{{unit}}.txt",
    benchmark:
        f"{BENCHDIR}/signal_qc/{{unit}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/signal_qc.py \
            --unit {wildcards.unit:q} \
            --plus {input.plus:q} \
            --minus {input.minus:q} \
            --gtf {input.gtf:q} \
            --chrom_sizes {input.sizes:q} \
            --peaks {input.peaks:q} \
            --output {output:q}

        sync {output:q}
        """


rule replicate_correlation:
    input:
        plus=lambda w: [f"{prefix(lib)}.+.bw" for lib in LIBRARIES_OF[w.experiment]],
        minus=lambda w: [f"{prefix(lib)}.-.bw" for lib in LIBRARIES_OF[w.experiment]],
        sizes=f"{REFDIR}/signal.chrom.sizes",
    output:
        f"{OUTDIR}/{{experiment}}/{{experiment}}.replicate_correlation.tsv",
    params:
        libraries=lambda w: LIBRARIES_OF[w.experiment],
    log:
        f"{LOGDIR}/replicate_correlation/{{experiment}}.txt",
    benchmark:
        f"{BENCHDIR}/replicate_correlation/{{experiment}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/replicate_correlation.py \
            --experiment {wildcards.experiment:q} \
            --libraries {params.libraries:q} \
            --plus {input.plus:q} \
            --minus {input.minus:q} \
            --chrom_sizes {input.sizes:q} \
            --output {output:q}

        sync {output:q}
        """


rule qc_table:
    input:
        samples=config["samples"],
        fastp=[f"{prefix(r)}.fastp.json" for r in RUNS],
        star=[f"{prefix(r)}.Log.final.out" for r in RUNS],
        umi_logs=[f"{prefix(lib)}.umi_tools.log" for lib in LIBRARIES],
        library_stats=[f"{prefix(lib)}.samtools_stats.txt" for lib in LIBRARIES],
        idxstats=[f"{prefix(lib)}.idxstats.tsv" for lib in LIBRARIES],
        signal_qc=[f"{prefix(u)}.signal_qc.tsv" for u in LIBRARIES + EXPERIMENTS],
        correlations=[f"{prefix(e)}.replicate_correlation.tsv" for e in REPLICATED],
        peaks=[bed for e in EXPERIMENTS for bed in peak_beds(e)],
    output:
        libraries=f"{OUTDIR}/qc/libraries.tsv",
        experiments=f"{OUTDIR}/qc/experiments.tsv",
    log:
        f"{LOGDIR}/qc_table.txt",
    benchmark:
        f"{BENCHDIR}/qc_table.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        python {workflow.basedir:q}/scripts/qc_table.py \
            --samples {input.samples:q} \
            --fastp {input.fastp:q} \
            --star {input.star:q} \
            --umi_logs {input.umi_logs:q} \
            --library_stats {input.library_stats:q} \
            --idxstats {input.idxstats:q} \
            --signal_qc {input.signal_qc:q} \
            --correlations {input.correlations:q} \
            --peaks {input.peaks:q} \
            --libraries {output.libraries:q} \
            --experiments {output.experiments:q}

        sync {output:q}
        """


rule multiqc:
    input:
        [f"{prefix(r)}.fastp.json" for r in RUNS],
        [f"{prefix(r)}.{s}_{m}_fastqc.zip" for r in RUNS for s in ["raw", "trimmed"] for m in ["R1", "R2"]],
        [f"{prefix(r)}.{f}" for r in RUNS for f in ["Log.final.out", "ReadsPerGene.out.tab"]],
        [f"{prefix(n)}.samtools_stats.txt" for n in RUNS + LIBRARIES],
        [f"{prefix(lib)}.umi_tools.log" for lib in LIBRARIES],
    output:
        html=f"{OUTDIR}/qc/multiqc.html",
        data=directory(f"{OUTDIR}/qc/multiqc_data"),
    params:
        outdir=f"{OUTDIR}/qc",
    log:
        f"{LOGDIR}/multiqc.txt",
    benchmark:
        f"{BENCHDIR}/multiqc.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        files=$(mktemp)
        trap 'rm -f "$files"' EXIT
        printf '%s\n' {input:q} > "$files"

        multiqc \
            --file-list "$files" \
            --outdir {params.outdir:q} \
            --filename multiqc.html \
            --force

        sync {output:q}
        """
