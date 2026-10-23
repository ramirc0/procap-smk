"""Join the per-run, per-library and per-experiment QC into two tables.

File names start with the unit: `<experiment>.<library>.<run>` for a run,
`<experiment>.<library>` for a library and `<experiment>` for an experiment.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import polars as pl

PEAK_TYPES = ["bidirectional", "divergent", "unidirectional"]


def build_parser():
    """Build the qc_table.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the qc_table.py flags.
    """
    p = argparse.ArgumentParser(description="Collect PRO-cap QC into two tables.")
    p.add_argument("--fastp", nargs="+", required=True, help="fastp JSON per run.")
    p.add_argument(
        "--star", nargs="+", required=True, help="STAR Log.final.out per run."
    )
    p.add_argument(
        "--umi_logs", nargs="+", required=True, help="umi_tools dedup log per library."
    )
    p.add_argument(
        "--library_stats",
        nargs="+",
        required=True,
        help="samtools stats per library BAM.",
    )
    p.add_argument(
        "--idxstats",
        nargs="+",
        required=True,
        help="samtools idxstats per library BAM.",
    )
    p.add_argument(
        "--signal_qc",
        nargs="+",
        required=True,
        help="signal_qc.py row per library and experiment.",
    )
    p.add_argument(
        "--correlations",
        nargs="*",
        default=[],
        help="replicate_correlation.py table per experiment.",
    )
    p.add_argument("--peaks", nargs="+", required=True, help="PINTS peak BEDs.")
    p.add_argument("--libraries", required=True, help="Library table to write.")
    p.add_argument("--experiments", required=True, help="Experiment table to write.")
    return p


def main(argv=None):
    """Parse every report and write the library and experiment tables.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    runs = pl.DataFrame([_run_row(f, s) for f, s in zip(args.fastp, args.star)])
    signal = pl.concat(
        [pl.read_csv(f, separator="\t") for f in args.signal_qc], how="diagonal_relaxed"
    )
    is_library = pl.col("unit").str.contains(".", literal=True)
    unit_parts = pl.col("unit").str.split_exact(".", 1)

    libraries = (
        runs.group_by("experiment", "library")
        .agg(
            pl.len().alias("runs"),
            pl.col(
                "raw_pairs", "trimmed_pairs", "star_input_pairs", "star_unique_pairs"
            ).sum(),
        )
        .join(
            pl.DataFrame(
                [
                    _library_row(*f)
                    for f in zip(args.umi_logs, args.library_stats, args.idxstats)
                ]
            ),
            on=["experiment", "library"],
        )
        .join(
            signal.filter(is_library).select(
                unit_parts.struct.field("field_0").alias("experiment"),
                unit_parts.struct.field("field_1").alias("library"),
                "gene_body_ratio",
            ),
            on=["experiment", "library"],
        )
        .with_columns(
            (100 * pl.col("star_unique_pairs") / pl.col("star_input_pairs")).alias(
                "star_unique_pct"
            ),
            (1 - pl.col("dedup_pairs") / pl.col("filtered_pairs")).alias(
                "duplication_rate"
            ),
        )
        .select(
            "experiment",
            "library",
            "runs",
            "raw_pairs",
            "trimmed_pairs",
            "star_input_pairs",
            "star_unique_pairs",
            "star_unique_pct",
            "filtered_pairs",
            "dedup_pairs",
            "duplication_rate",
            "signal_pairs",
            "spikein_pairs",
            "rdna_pairs",
            "chrM_pairs",
            "median_rna_length",
            "gene_body_ratio",
        )
        .sort("experiment", "library")
    )
    libraries.write_csv(args.libraries, separator="\t")

    experiments = (
        libraries.group_by("experiment")
        .agg(pl.len().alias("libraries"), pl.col("dedup_pairs").sum())
        .join(
            signal.filter(~is_library).select(
                pl.col("unit").alias("experiment"),
                (pl.col("plus_signal") + pl.col("minus_signal")).alias("signal_pairs"),
                "gene_body_ratio",
                "fraction_in_peaks",
            ),
            on="experiment",
        )
        .join(_peak_counts(args.peaks), on="experiment")
        .join(_correlations(args.correlations), on="experiment", how="left")
        .sort("experiment")
    )
    experiments.write_csv(args.experiments, separator="\t")


def _unit(path, suffix):
    """Return the run, library or experiment a QC file belongs to.

    Parameters
    ----------
    path : str
        QC file path, `.../<rule>/<unit><suffix>`.
    suffix : str
        File name ending after the unit.

    Returns
    -------
    list of str
        The unit split into experiment, then library and run if present.
    """
    return Path(path).name.removesuffix(suffix).split(".")


def _run_row(fastp, star):
    """Read the pair counts of one run from fastp and STAR.

    Parameters
    ----------
    fastp : str
        The fastp JSON report.
    star : str
        The STAR Log.final.out.

    Returns
    -------
    dict
        Experiment, library and pair counts.
    """
    summary = json.loads(Path(fastp).read_text())["summary"]
    pairs = (
        line.split("|", 1)
        for line in Path(star).read_text().splitlines()
        if "|" in line
    )
    log = {k.strip(): v.strip() for k, v in pairs}
    experiment, library, _ = _unit(fastp, ".json")
    return {
        "experiment": experiment,
        "library": library,
        "raw_pairs": summary["before_filtering"]["total_reads"] // 2,
        "trimmed_pairs": summary["after_filtering"]["total_reads"] // 2,
        "star_input_pairs": int(log["Number of input reads"]),
        "star_unique_pairs": int(log["Uniquely mapped reads number"]),
    }


def _library_row(umi_log, stats, idxstats):
    """Read the pair counts of one library from umi_tools and samtools.

    Parameters
    ----------
    umi_log : str
        The umi_tools dedup log.
    stats : str
        Output of samtools stats on the deduplicated BAM.
    idxstats : str
        Output of samtools idxstats on the deduplicated BAM.

    Returns
    -------
    dict
        Experiment, library, pair counts and median RNA length.
    """
    text = Path(umi_log).read_text()
    filtered = int(re.search(r"Input Reads: (\d+)", text).group(1))
    summary = {}
    sizes, pairs = [], []
    for line in Path(stats).read_text().splitlines():
        fields = line.split("\t")
        if fields[0] == "SN":
            summary[fields[1].rstrip(":")] = fields[2]
        elif fields[0] == "IS":
            sizes.append(int(fields[1]))
            pairs.append(int(fields[2]))
    contigs = pl.read_csv(
        idxstats,
        separator="\t",
        has_header=False,
        new_columns=["contig", "length", "mapped", "unmapped"],
    )
    contig = pl.col("contig")
    experiment, library = _unit(stats, ".txt")
    return {
        "experiment": experiment,
        "library": library,
        "filtered_pairs": filtered,
        "dedup_pairs": int(summary["raw total sequences"]) // 2,
        "signal_pairs": _pairs_on(contigs, contig.str.contains(r"^chr([0-9]+|X|Y)$")),
        "spikein_pairs": _pairs_on(contigs, contig.str.starts_with("dm6_")),
        "rdna_pairs": _pairs_on(contigs, contig == "U13369.1"),
        "chrM_pairs": _pairs_on(contigs, contig == "chrM"),
        "median_rna_length": _weighted_median(np.array(sizes), np.array(pairs)),
    }


def _pairs_on(contigs, selected):
    """Count the read pairs mapped to some contigs.

    Parameters
    ----------
    contigs : polars.DataFrame
        The samtools idxstats table.
    selected : polars.Expr
        Filter that picks the contigs.

    Returns
    -------
    int
        Mapped reads on those contigs, halved.
    """
    return int(contigs.filter(selected)["mapped"].sum()) // 2


def _weighted_median(values, weights):
    """Return the median of values repeated by weights.

    Parameters
    ----------
    values : numpy.ndarray
        Sorted values.
    weights : numpy.ndarray
        Count of each value.

    Returns
    -------
    float or None
        The median, or None without weight.
    """
    if weights.sum() == 0:
        return None
    cumulative = np.cumsum(weights)
    return float(values[np.searchsorted(cumulative, cumulative[-1] / 2)])


def _peak_counts(paths):
    """Count the peaks of each type per experiment.

    Parameters
    ----------
    paths : list of str
        PINTS BEDs, `<experiment>.<type>.bed`.

    Returns
    -------
    polars.DataFrame
        One row per experiment, one column per peak type.
    """
    rows = {}
    for path in paths:
        kind = next(t for t in PEAK_TYPES if Path(path).name.endswith(f".{t}.bed"))
        (experiment,) = _unit(path, f".{kind}.bed")
        with open(path) as f:
            rows.setdefault(experiment, {})[f"{kind}_peaks"] = sum(1 for _ in f)
    return pl.DataFrame([{"experiment": e, **counts} for e, counts in rows.items()])


def _correlations(paths):
    """Summarize the library correlations of each experiment.

    Parameters
    ----------
    paths : list of str
        Tables written by replicate_correlation.py.

    Returns
    -------
    polars.DataFrame
        Columns experiment, min_replicate_r, median_replicate_r.
    """
    schema = {
        "experiment": pl.String,
        "min_replicate_r": pl.Float64,
        "median_replicate_r": pl.Float64,
    }
    if not paths:
        return pl.DataFrame(schema=schema)
    return (
        pl.concat([pl.read_csv(p, separator="\t") for p in paths])
        .group_by("experiment")
        .agg(
            pl.col("pearson_r").min().alias("min_replicate_r"),
            pl.col("pearson_r").median().alias("median_replicate_r"),
        )
    )


if __name__ == "__main__":
    main()
