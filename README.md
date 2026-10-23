# procap-smk

A Snakemake port of ENCODE's PRO-cap pipeline v2, "Yu Lab Nuclear Run-On Assay
Pipeline v2" (ENCPL998WCJ), run with current tool versions. It goes from FASTQ
to filtered BAMs, RNA 5' end bigWigs, PINTS peaks and QC tables.

## Steps

| Unit | Step | Tool |
|---|---|---|
| run (FASTQ pair) | UMI extraction (6 nt per read) and trimming | fastp |
| run | alignment to GRCh38 + rDNA `U13369.1` (+ dm6 for spike-in runs) | STAR |
| library | merge runs; keep MAPQ 255 proper pairs with no spliced mate | samtools |
| library | UMI deduplication | umi_tools |
| library, experiment | read 2 5' ends per strand, one count per pair | figwig |
| experiment | bidirectional, divergent and unidirectional peaks | PINTS |

The flags are ENCODE v2's, recovered from the rmsp-library template
`PROcapAnalysis/20240501` and the `@PG` lines of ENCODE's BAMs. The fastp
adapters are the exception: ENCODE's trimmed reads end at the Illumina RA3 and
RA5 adapters, not at the template's defaults. The bigWigs hold raw counts on
chr1-22, X and Y. Minus-strand values are positive.

## Layout

```
.tests/integration/      config + toy inputs for the end-to-end test; make_data.sh rebuilds them
.tests/encode/           compare_encode.py: agreement with ENCODE's released files
config/
  config.yaml.template   references and parameters; copy to config.yaml (gitignored)
  samples.example.tsv    one row per FASTQ pair: run, library, experiment, r1, r2, spikein
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, run sheet, unit maps, helpers
    refs.smk             combined genome FASTAs, STAR indexes, signal chrom sizes
    align.smk            fastp, STAR
    filter.smk           unique-pair filter, UMI deduplication
    signal.smk           bigWigs
    peaks.smk            PINTS
    qc.smk               FastQC, samtools stats, signal QC, QC tables, MultiQC
  scripts/               argparse scripts, each with build_parser() + main()
  envs/pixi.toml         rule env (pixi workspace)
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + test and lint tasks
pyproject.toml           ruff + numpydoc config
```

## Inputs

- **Run sheet.** One row per FASTQ pair. Runs of one `library` are merged
  before deduplication. A library label such as `1_2` is unique within its
  `experiment`. Libraries of one experiment are pooled for peak calling.
  `spikein` is `dm6` for runs aligned to the Drosophila spike-in index, else
  empty. IDs must not contain dots.
- **FASTQs** must be on disk. The workflow downloads nothing.
- **References** are local paths in the config: the GRCh38 no-alt analysis set,
  the rDNA repeat `U13369.1`, UCSC dm6 and GENCODE v50. FASTAs may be gzipped;
  the GTF must be plain text.

## Outputs

Under `results/<run_id>/<experiment>/<rule>/<unit>`. A run's unit is
`<experiment>.<library>.<run>`, a library's `<experiment>.<library>`, an
experiment's `<experiment>`. Logs and benchmarks follow the same layout.

| Rule | Unit | Files |
|---|---|---|
| `fastp` | run | `.json`, `.html` |
| `fastqc` | run | `.{raw,trimmed}_R{1,2}_fastqc.{zip,html}` |
| `star_align` | run | `.bam`, `.Log.final.out`, `.ReadsPerGene.out.tab` |
| `dedup` | library | `.bam`, `.bam.bai`, `.umi_tools.log` |
| `bigwig` | library, experiment | `.+.bw`, `.-.bw` |
| `pints` | experiment | `.{bidirectional,divergent,unidirectional}.bed` |
| `samtools_stats` | run, library | `.txt` |
| `idxstats` | library | `.tsv` |
| `signal_qc` | library, experiment | `.tsv` |
| `replicate_correlation` | experiment | `.tsv` |

`results/<run_id>/qc/` holds `libraries.tsv`, `experiments.tsv` and the MultiQC
report. The gene-body ratio follows the paper's description (Shah et al.,
bioRxiv 10.1101/2025.09.24.676871). Replicate correlation is Pearson r of
log10 counts in 50 bp bins that are nonzero in both libraries.

## Running

```bash
cp config/config.yaml.template config/config.yaml  # then edit
pixi shell                                         # launcher env
snakemake --profile profiles/slurm                 # SLURM
snakemake -n -p --profile profiles/local           # dry run
pixi run test                                      # end-to-end test
```

## Environment

The launcher env (`pixi.toml` at the root) has Snakemake main, the pixi
software-deployment plugin and the SLURM executor plugins. The rule env is a
second pixi workspace, `workflow/envs/`, and its `pixi.lock` pins every package.
An `onstart` hook installs it before the first job. Any change to the manifest
or lock reruns every rule.
