"""Compare procap-smk outputs with ENCODE's released PRO-cap v2 files.

Per run: the unfiltered BAM's @SQ lines and record count. Per library: the
filtered BAM's pair count, and per-base and 50 bp bin agreement of each strand's
bigWig (ENCODE's minus strand is negative, ours positive). Per experiment: the
count and reciprocal overlap of each peak type. ENCODE BAMs stream over HTTPS;
bigWigs and peaks come from the local ENCODE mirror.
"""

import argparse
import gzip
import subprocess
from pathlib import Path

import numpy as np
import polars as pl
from figwig import BigWigReader

PORTAL = "https://www.encodeproject.org/files"
PEAK_TYPES = ["bidirectional", "divergent", "unidirectional"]


def build_parser():
    """Build the compare_encode.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the compare_encode.py flags.
    """
    p = argparse.ArgumentParser(description="Compare procap-smk outputs with ENCODE.")
    p.add_argument("--metadata", required=True, help="ENCODE metadata Parquet dir.")
    p.add_argument("--samples", required=True, help="Run sheet of the validation run.")
    p.add_argument("--results", required=True, help="Output dir, results/<run_id>.")
    p.add_argument("--encode_root", required=True, help="Local ENCODE mirror.")
    p.add_argument("--chrom_sizes", required=True, help="Signal chrom sizes.")
    p.add_argument("-o", "--output", required=True, help="Long-format TSV to write.")
    return p


def main(argv=None):
    """Compare every run, library and experiment of the sheet.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    sheet = pl.read_csv(args.samples, separator="\t", infer_schema=False)
    sizes = pl.read_csv(args.chrom_sizes, separator="\t", has_header=False).rows()
    files = (
        pl.scan_parquet(f"{args.metadata}/file.parquet")
        .filter(
            pl.col("dataset").is_in(
                [f"/experiments/{e}/" for e in sheet["experiment"].unique()]
            ),
            pl.col("status") == "released",
            pl.col("analysis_step_version").str.contains("v2-step"),
        )
        .select(
            pl.col("dataset").str.extract(r"(ENCSR\w+)").alias("experiment"),
            "accession",
            "output_type",
            "technical_replicates",
            "derived_from",
        )
        .collect()
    )
    out = Path(args.results)
    rows = []
    for run, library, experiment, r1 in sheet.select(
        "run", "library", "experiment", "r1"
    ).iter_rows():
        acc = Path(r1).name.removesuffix(".fastq.gz")
        encode = files.filter(
            pl.col("output_type") == "unfiltered alignments",
            pl.col("derived_from").str.contains(f"/files/{acc}/"),
        )["accession"].item()
        ours = out / experiment / "star_align" / f"{experiment}.{library}.{run}.bam"
        rows.append(
            (
                "run",
                run,
                "sq_identical",
                _sq(str(ours)) == _sq(_url(encode, "bam")),
                True,
            )
        )
        rows.append(
            ("run", run, "records", _count(str(ours)), _count(_url(encode, "bam")))
        )
    for library, experiment in (
        sheet.select("library", "experiment").unique().iter_rows()
    ):
        mine = files.filter(
            pl.col("experiment") == experiment,
            pl.col("technical_replicates") == library,
        )
        stem = f"{experiment}.{library}"
        encode_bam = mine.filter(pl.col("output_type") == "alignments")[
            "accession"
        ].item()
        rows.append(
            (
                "library",
                stem,
                "pairs",
                _count(str(out / experiment / "dedup" / f"{stem}.bam")) // 2,
                _count(_url(encode_bam, "bam")) // 2,
            )
        )
        for strand, kind in [("+", "plus"), ("-", "minus")]:
            encode_bw = mine.filter(pl.col("output_type").str.starts_with(kind))[
                "accession"
            ].item()
            for metric, ours, theirs in _bigwig_agreement(
                out / experiment / "bigwig" / f"{stem}.{strand}.bw",
                Path(args.encode_root, experiment, f"{encode_bw}.bigWig"),
                sizes,
            ):
                rows.append(("library", stem, f"{kind}_{metric}", ours, theirs))
    for experiment in sheet["experiment"].unique().sort():
        for kind in PEAK_TYPES:
            encode_bed = files.filter(
                pl.col("experiment") == experiment,
                pl.col("output_type") == f"{kind} peaks",
            )["accession"].item()
            ours = out / experiment / "pints" / f"{experiment}.{kind}.bed"
            theirs = Path(args.encode_root, experiment, f"{encode_bed}.bed.gz")
            rows += [
                (
                    "experiment",
                    experiment,
                    f"{kind}_peaks",
                    _lines(ours),
                    _lines(theirs),
                ),
                (
                    "experiment",
                    experiment,
                    f"{kind}_overlapped",
                    _overlapped(ours, theirs),
                    _overlapped(theirs, ours),
                ),
            ]
    report = pl.DataFrame(
        [tuple(str(v) for v in r) for r in rows],
        schema=["level", "unit", "metric", "ours", "encode"],
        orient="row",
    )
    report.write_csv(args.output, separator="\t")
    with pl.Config(tbl_rows=-1, tbl_width_chars=200):
        print(report)


def _url(accession, extension):
    """Return the ENCODE portal download URL of a file.

    Parameters
    ----------
    accession : str
        File accession.
    extension : str
        File extension, such as `bam`.

    Returns
    -------
    str
        The download URL.
    """
    return f"{PORTAL}/{accession}/@@download/{accession}.{extension}"


def _sq(bam):
    """Return the @SQ lines of a BAM header.

    Parameters
    ----------
    bam : str
        Path or URL.

    Returns
    -------
    list of str
        The @SQ lines, in order.
    """
    header = subprocess.run(
        ["samtools", "view", "-H", bam], check=True, capture_output=True, text=True
    ).stdout
    return [line for line in header.splitlines() if line.startswith("@SQ")]


def _count(bam):
    """Count the records of a BAM.

    Parameters
    ----------
    bam : str
        Path or URL.

    Returns
    -------
    int
        Number of records.
    """
    result = subprocess.run(
        ["samtools", "view", "-c", bam], check=True, capture_output=True, text=True
    )
    return int(result.stdout)


def _lines(path):
    """Count the lines of a text file, gzipped or not.

    Parameters
    ----------
    path : pathlib.Path
        File path.

    Returns
    -------
    int
        Number of lines.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as f:
        return sum(1 for _ in f)


def _overlapped(a, b):
    """Return the fraction of intervals in `a` that overlap any interval in `b`.

    Parameters
    ----------
    a, b : pathlib.Path
        BED files, gzipped or not.

    Returns
    -------
    float
        Overlapping fraction of `a`.
    """
    hits = subprocess.run(
        ["bedtools", "intersect", "-u", "-a", str(a), "-b", str(b)],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return len(hits.splitlines()) / max(_lines(a), 1)


def _bigwig_agreement(ours, theirs, sizes):
    """Compare two tracks of one strand base by base and in 50 bp bins.

    Parameters
    ----------
    ours : pathlib.Path
        Our bigWig.
    theirs : pathlib.Path
        ENCODE's bigWig; its values are compared as absolute values.
    sizes : list of tuple of (str, int)
        Chromosomes and lengths.

    Returns
    -------
    list of tuple
        Metric name, our value, ENCODE's value.
    """
    a_reader, b_reader = BigWigReader(str(ours)), BigWigReader(str(theirs))
    total_a = total_b = 0.0
    covered = identical = 0
    bins_a, bins_b = [], []
    for chrom, length in sizes:
        a = _read(a_reader, chrom, length)
        b = np.abs(_read(b_reader, chrom, length))
        total_a += a.sum(dtype=np.float64)
        total_b += b.sum(dtype=np.float64)
        either = (a != 0) | (b != 0)
        covered += int(either.sum())
        identical += int((a[either] == b[either]).sum())
        n = length // 50
        bins_a.append(a[: n * 50].reshape(n, 50).sum(axis=1, dtype=np.float64))
        bins_b.append(b[: n * 50].reshape(n, 50).sum(axis=1, dtype=np.float64))
    r = np.corrcoef(np.concatenate(bins_a), np.concatenate(bins_b))[0, 1]
    return [
        ("total", total_a, total_b),
        ("identical_bases", identical / max(covered, 1), 1.0),
        ("bin50_pearson", r, 1.0),
    ]


def _read(reader, chrom, length):
    """Read one chromosome of a track, zeros where it has no data.

    Parameters
    ----------
    reader : figwig.BigWigReader
        Open track.
    chrom : str
        Chromosome name.
    length : int
        Chromosome length.

    Returns
    -------
    numpy.ndarray
        Signal per base, float32.
    """
    if chrom not in reader.chrom_sizes:
        return np.zeros(length, dtype=np.float32)
    return reader.read(chrom, [0], length)[0]


if __name__ == "__main__":
    main()
