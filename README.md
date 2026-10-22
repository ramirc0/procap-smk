# procap-smk

A Snakemake port of ENCODE's PRO-cap pipeline v2, "Yu Lab Nuclear Run-On Assay
Pipeline v2" (ENCPL998WCJ), run with current tool versions.

## Layout

```
.tests/integration/      config + toy inputs for the end-to-end test
config/
  config.yaml.template   references and parameters; copy to config.yaml (gitignored)
  samples.example.tsv    one row per FASTQ pair: run, library, experiment, r1, r2, spikein
workflow/
  Snakefile              includes the rules + rule all
  rules/
    common.smk           config, run sheet, unit maps, helpers
    refs.smk             combined genome FASTAs, STAR indexes, signal chrom sizes
  envs/pixi.toml         rule env (pixi workspace)
profiles/
  local/config.yaml      local execution
  slurm/config.yaml      SLURM (per-rule resources, retried with more memory)
pixi.toml                launcher env + test and lint tasks
pyproject.toml           ruff + numpydoc config
```

## Inputs

- **Run sheet.** One row per FASTQ pair. Runs of one `library` are merged
  before deduplication. Libraries of one `experiment` are pooled for peak
  calling. `spikein` is `dm6` for runs aligned to the Drosophila spike-in
  index, else empty.
- **FASTQs** must be on disk. The workflow downloads nothing.
- **References** are local paths in the config. The genome is the GRCh38 no-alt
  analysis set, the rDNA repeat is `U13369.1`, the spike-in is UCSC dm6, and the
  annotation is GENCODE v50.

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
