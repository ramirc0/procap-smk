#!/usr/bin/env bash
# Build the toy inputs of `pixi run test` from ENCODE data. Run once from the
# repository root, in the rule env:
#   pixi run --manifest-path workflow/envs/pixi.toml bash .tests/integration/make_data.sh
# It reads the procap-atlas FASTQs and streams ENCODE's unfiltered BAMs to pick
# the read pairs ENCODE aligned to the toy genome.
set -euo pipefail

static=/data/zusers/ramirezc/static
fastqs=/home/ramirezc/Projects/jackfruit/data/procap-atlas/encode/fastqs
out=.tests/integration/resources
start=6000000   # chr12 slice [start, end) around GAPDH, renamed chr12
end=7000000
max_pairs=15000 # per library in the slice; rDNA and dm6 take a tenth each
mkdir -p "$out/reads"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

samtools faidx "$static/hg38/genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta" \
    "chr12:$((start + 1))-$end" | sed '1s/.*/>chr12/' > "$out/genome.fa"
install -m 644 "$static/hg38/genome/U13369.1.fasta" "$out/rdna.fa"
samtools faidx "$static/dm6/genome/dm6.fa" chr4:1-200000 | sed '1s/.*/>chr4/' > "$out/spikein.fa"

# GENCODE v50 genes that lie inside the slice, shifted to slice coordinates.
# The workflow reads only gene, transcript and exon lines.
gtf="$static/hg38/annotation/gencode.v50.annotation.gtf"
awk -F '\t' -v s="$start" -v e="$end" '$1 == "chr12" && $3 == "gene" && $4 > s && $5 <= e {
    match($9, /gene_id "[^"]+"/); print substr($9, RSTART, RLENGTH) }' "$gtf" > "$tmp/genes"
awk -F '\t' -v OFS='\t' -v s="$start" 'NR == FNR { want[$0]; next }
    $1 == "chr12" && $3 ~ /^(gene|transcript|exon)$/ {
        match($9, /gene_id "[^"]+"/); if (substr($9, RSTART, RLENGTH) in want) { $4 -= s; $5 -= s; print } }' \
    "$tmp/genes" "$gtf" > "$out/annotation.gtf"

# FASTQ names of the pairs an ENCODE BAM placed in the toy genome. fastp
# appended ":<UMI1>_<UMI2>" to them, so that suffix comes off.
names() {
    samtools view -F 0x904 "https://www.encodeproject.org/files/$1/@@download/$1.bam" \
        | awk -v s="$start" -v e="$end" -v n="$max_pairs" '
            { sub(/:[ACGTN]+_[ACGTN]+$/, "", $1) }
            $1 in seen { next }
            $3 == "chr12" && $4 > s && $4 <= e - 500 && human < n { seen[$1]; human++; print $1; next }
            $3 == "U13369.1" && rdna < n / 10 { seen[$1]; rdna++; print $1; next }
            $3 == "dm6_chr4" && $4 <= 199500 && fly < n / 10 { seen[$1]; fly++; print $1 }'
}

# Write the reads of the listed pairs from one FASTQ, gzipped.
extract() {
    zcat "$2" | awk 'NR == FNR { want[$1]; next }
        FNR % 4 == 1 { split($1, a, " "); keep = (substr(a[1], 2) in want) }
        keep' "$1" - | gzip -n > "$3"
}

# toyA: spinal cord ENCSR353JON (no spike-in), split into two runs of one library.
names ENCFF348WGN > "$tmp/a"
awk 'NR % 2' "$tmp/a" > "$tmp/a1"
awk '!(NR % 2)' "$tmp/a" > "$tmp/a2"
for run in 1 2; do
    extract "$tmp/a$run" "$fastqs/ENCFF621YNC.fastq.gz" "$out/reads/toyA_r${run}_R1.fastq.gz"
    extract "$tmp/a$run" "$fastqs/ENCFF563BPL.fastq.gz" "$out/reads/toyA_r${run}_R2.fastq.gz"
done

# toyB: HCT116 ENCSR197JLW (dm6 spike-in), one run from each of two libraries.
names ENCFF845FEM > "$tmp/b1"
names ENCFF258GQP > "$tmp/b2"
extract "$tmp/b1" "$fastqs/ENCFF254NZI.fastq.gz" "$out/reads/toyB_r1_R1.fastq.gz"
extract "$tmp/b1" "$fastqs/ENCFF982HII.fastq.gz" "$out/reads/toyB_r1_R2.fastq.gz"
extract "$tmp/b2" "$fastqs/ENCFF546OUD.fastq.gz" "$out/reads/toyB_r2_R1.fastq.gz"
extract "$tmp/b2" "$fastqs/ENCFF639FYO.fastq.gz" "$out/reads/toyB_r2_R2.fastq.gz"

wc -l "$tmp"/a "$tmp"/b1 "$tmp"/b2
