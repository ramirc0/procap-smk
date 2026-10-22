# Combined genome FASTAs, STAR indexes and the signal chrom sizes. ENCODE aligned
# spike-in runs to hg38_rDNA_dm6 (dm6 contigs prefixed dm6_), the rest to hg38_rDNA.


# Run-independent: built once under results/refs/, reused across runs.
rule genome_fasta:
    input:
        genome=REFERENCES["genome"],
        rdna=REFERENCES["rdna"],
        spikein=REFERENCES["spikein"],
    output:
        temp(f"{REFDIR}/{{genome}}/genome.fa"),
    params:
        spikein=lambda w, input: input.spikein if w.genome == "hg38_rDNA_dm6" else "",
    log:
        f"{LOGDIR}/genome_fasta/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/genome_fasta/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        {{
            zcat -f {input.genome:q} {input.rdna:q}
            if [[ -n {params.spikein:q} ]]; then
                zcat -f {params.spikein:q} | sed 's/^>/>dm6_/'
            fi
        }} > {output:q}

        sync {output:q}
        """


rule star_index:
    input:
        fasta=f"{REFDIR}/{{genome}}/genome.fa",
        gtf=REFERENCES["gtf"],
    output:
        directory(f"{REFDIR}/{{genome}}/star"),
    params:
        sa_index_nbases=lambda w: config["star"]["sa_index_nbases"],
        ram_bytes=lambda w, resources: resources.mem_mb * 1024**2,
    log:
        f"{LOGDIR}/star_index/{{genome}}.txt",
    benchmark:
        f"{BENCHDIR}/star_index/{{genome}}.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        ulimit -n 65536
        mkdir -p {output:q}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads:q} \
            --genomeDir {output:q} \
            --genomeFastaFiles {input.fasta:q} \
            --sjdbGTFfile {input.gtf:q} \
            --genomeSAindexNbases {params.sa_index_nbases:q} \
            --limitGenomeGenerateRAM {params.ram_bytes:q} \
            --outFileNamePrefix {output:q}/

        sync {output:q}
        """


# chr1-22, X and Y, in FASTA order: the contigs ENCODE's signal tracks keep.
rule chrom_sizes:
    input:
        REFERENCES["genome"],
    output:
        f"{REFDIR}/signal.chrom.sizes",
    log:
        f"{LOGDIR}/chrom_sizes.txt",
    benchmark:
        f"{BENCHDIR}/chrom_sizes.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        zcat -f {input:q} \
            | awk '/^>/ {{ if (keep) print name "\t" len; name = substr($1, 2); len = 0; keep = name ~ /^chr([0-9]+|X|Y)$/; next }}
                   keep {{ len += length($0) }}
                   END {{ if (keep) print name "\t" len }}' \
            > {output:q}

        sync {output:q}
        """
