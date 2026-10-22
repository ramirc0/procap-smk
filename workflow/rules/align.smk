# Per run: UMI extraction and trimming, then alignment. Flags are ENCODE v2's
# (rmsp-library PROcapAnalysis 20240501 and the STAR @PG line of its BAMs).


rule fastp:
    input:
        r1=lambda w: R1_OF[w.run],
        r2=lambda w: R2_OF[w.run],
    output:
        r1=temp(f"{OUTDIR}/{{run}}/{{run}}_R1.fastq.gz"),
        r2=temp(f"{OUTDIR}/{{run}}/{{run}}_R2.fastq.gz"),
        html=f"{OUTDIR}/{{run}}/{{run}}.fastp.html",
        json=f"{OUTDIR}/{{run}}/{{run}}.fastp.json",
    log:
        f"{LOGDIR}/fastp/{{run}}.txt",
    benchmark:
        f"{BENCHDIR}/fastp/{{run}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        fastp \
            --in1 {input.r1:q} \
            --in2 {input.r2:q} \
            --out1 {output.r1:q} \
            --out2 {output.r2:q} \
            --html {output.html:q} \
            --json {output.json:q} \
            --thread {threads:q} \
            --umi \
            --umi_loc per_read \
            --umi_len 6 \
            --adapter_sequence ACGGAT \
            --adapter_sequence_r2 CGGCTT \
            --overlap_len_require 18 \
            --length_required 18 \
            --low_complexity_filter \
            --trim_poly_g \
            --overrepresentation_analysis \
            --correction

        sync {output:q}
        """


rule star_align:
    input:
        index=lambda w: f"{REFDIR}/{genome_of(w.run)}/star",
        r1=f"{OUTDIR}/{{run}}/{{run}}_R1.fastq.gz",
        r2=f"{OUTDIR}/{{run}}/{{run}}_R2.fastq.gz",
    output:
        bam=f"{OUTDIR}/{{run}}/{{run}}.unfiltered.bam",
        log_final=f"{OUTDIR}/{{run}}/{{run}}.Log.final.out",
        gene_counts=f"{OUTDIR}/{{run}}/{{run}}.ReadsPerGene.out.tab",
    log:
        f"{LOGDIR}/star_align/{{run}}.txt",
    benchmark:
        f"{BENCHDIR}/star_align/{{run}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536
        out=$(mktemp -d)
        trap 'rm -rf "$out"' EXIT

        STAR \
            --runThreadN {threads:q} \
            --genomeDir {input.index:q} \
            --readFilesIn {input.r1:q} {input.r2:q} \
            --readFilesCommand zcat \
            --outFileNamePrefix "$out/" \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --outFilterMultimapNmax 10 \
            --outFilterMultimapScoreRange 0 \
            --outFilterMismatchNmax 1 \
            --alignIntronMax 1000 \
            --alignMatesGapMax 1000 \
            --quantMode GeneCounts

        mv "$out/Aligned.sortedByCoord.out.bam" {output.bam:q}
        mv "$out/Log.final.out" {output.log_final:q}
        mv "$out/ReadsPerGene.out.tab" {output.gene_counts:q}

        sync {output:q}
        """
