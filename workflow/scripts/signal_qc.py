"""Signal QC for one PRO-cap track pair: totals, gene-body ratio, signal in peaks.

The gene-body ratio follows Shah et al. (bioRxiv 10.1101/2025.09.24.676871),
Methods "Data preprocessing". Genes are the protein-coding transcripts of at
least 2 kb whose gene overlaps no other gene. The top tenth by signal in
[TSS-100, TSS+1000) is kept, scoring each gene by its best transcript. Over
those genes, TSS regions are the first 500 bp of each transcript and TTS
regions the last 500 bp. The gene body is every transcript base outside both.
The ratio is (gene-body signal / length) / (TSS signal / length + gene-body
signal / length). The biodatatools code behind ENCODE's values ranks genes on
whole transcripts instead, so these values differ from theirs.
"""

import argparse

import numpy as np
import polars as pl
from figwig import BigWigReader

TSS_UPSTREAM = 100
TSS_DOWNSTREAM = 1000
END_REGION = 500
MIN_TRANSCRIPT = 2000
TOP_FRACTION = 0.1


def build_parser():
    """Build the signal_qc.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the signal_qc.py flags.
    """
    p = argparse.ArgumentParser(description="Signal QC of a PRO-cap track pair.")
    p.add_argument("--unit", required=True, help="Library or experiment ID.")
    p.add_argument("--plus", required=True, help="Plus-strand bigWig.")
    p.add_argument(
        "--minus", required=True, help="Minus-strand bigWig (positive values)."
    )
    p.add_argument("--gtf", required=True, help="GENCODE GTF, plain text.")
    p.add_argument("--chrom_sizes", required=True, help="Chromosomes to measure.")
    p.add_argument(
        "--peaks",
        nargs="*",
        default=[],
        help="Peak BEDs; adds the fraction of signal in peaks.",
    )
    p.add_argument("-o", "--output", required=True, help="One-row TSV to write.")
    return p


def main(argv=None):
    """Read the tracks chromosome by chromosome and write the QC row.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    sizes = dict(
        pl.read_csv(args.chrom_sizes, separator="\t", has_header=False).iter_rows()
    )
    readers = {"+": BigWigReader(args.plus), "-": BigWigReader(args.minus)}
    transcripts = _qc_transcripts(args.gtf, sizes)
    peaks = _peaks(args.peaks, sizes)

    totals = {"+": 0.0, "-": 0.0}
    scores = []
    for chrom, length in sizes.items():
        for strand, reader in readers.items():
            signal = _chromosome(reader, chrom, length)
            totals[strand] += float(signal.sum(dtype=np.float64))
            on = transcripts.filter(
                pl.col("chrom") == chrom, pl.col("strand") == strand
            )
            scores.append(on.with_columns(score=_window_sums(signal, on)))
    scores = pl.concat(scores)
    top = (
        scores.group_by("gene_id")
        .agg(pl.col("score").max())
        .sort("score", "gene_id", descending=[True, False])
        .head(max(1, int(TOP_FRACTION * scores["gene_id"].n_unique())))
    )
    chosen = transcripts.join(top.select("gene_id"), on="gene_id")

    sums = np.zeros(4)
    in_peaks = 0.0
    for chrom, length in sizes.items():
        for strand, reader in readers.items():
            signal = _chromosome(reader, chrom, length)
            on = chosen.filter(pl.col("chrom") == chrom, pl.col("strand") == strand)
            sums += _region_sums(signal, on)
            if peaks is not None:
                mask = _mask(length, peaks.filter(pl.col("chrom") == chrom))
                in_peaks += float(signal[mask].sum(dtype=np.float64))
    tss, tss_len, body, body_len = sums
    density = body / body_len if body_len else np.nan
    total = totals["+"] + totals["-"]
    pl.DataFrame(
        {
            "unit": [args.unit],
            "plus_signal": [totals["+"]],
            "minus_signal": [totals["-"]],
            "ratio_genes": [top.height],
            "gene_body_ratio": [
                density / (tss / tss_len + density) if tss_len else np.nan
            ],
            "fraction_in_peaks": [
                in_peaks / total if peaks is not None and total else None
            ],
        }
    ).write_csv(args.output, separator="\t")


def _qc_transcripts(gtf, sizes):
    """Return the protein-coding transcripts of genes that overlap no other gene.

    Parameters
    ----------
    gtf : str
        GENCODE GTF path.
    sizes : dict of str to int
        Chromosomes to keep.

    Returns
    -------
    polars.DataFrame
        Columns chrom, start, end (0-based, half-open), strand, gene_id.
    """
    rows = (
        pl.scan_csv(
            gtf,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            quote_char=None,
            schema={
                "chrom": pl.String,
                "source": pl.String,
                "feature": pl.String,
                "start": pl.Int64,
                "end": pl.Int64,
                "score": pl.String,
                "strand": pl.String,
                "frame": pl.String,
                "attributes": pl.String,
            },
        )
        .filter(
            pl.col("feature").is_in(["gene", "transcript"]),
            pl.col("chrom").is_in(list(sizes)),
        )
        .select(
            "chrom",
            (pl.col("start") - 1).alias("start"),
            "end",
            "strand",
            "feature",
            pl.col("attributes").str.extract(r'gene_id "([^"]+)"').alias("gene_id"),
            pl.col("attributes")
            .str.extract(r'transcript_type "([^"]+)"')
            .alias("type"),
        )
        .collect()
    )
    genes = rows.filter(pl.col("feature") == "gene").sort("chrom", "start")
    # A gene overlaps another when it starts before the furthest end seen so far
    # on its chromosome, or ends after the next gene starts.
    reach = pl.col("end").cum_max().shift(1).over("chrom")
    following = pl.col("start").shift(-1).over("chrom")
    isolated = genes.filter(
        (pl.col("start") >= reach).fill_null(True)
        & (pl.col("end") <= following).fill_null(True)
    )
    return rows.filter(
        pl.col("feature") == "transcript",
        pl.col("type") == "protein_coding",
        pl.col("end") - pl.col("start") >= MIN_TRANSCRIPT,
        pl.col("gene_id").is_in(isolated["gene_id"].implode()),
    ).select("chrom", "start", "end", "strand", "gene_id")


def _peaks(paths, sizes):
    """Return the peak intervals of every BED, or None without peaks.

    Parameters
    ----------
    paths : list of str
        Peak BED paths.
    sizes : dict of str to int
        Chromosomes to keep.

    Returns
    -------
    polars.DataFrame or None
        Columns chrom, start, end.
    """
    if not paths:
        return None
    frames = [
        pl.read_csv(p, separator="\t", has_header=False, infer_schema=False).select(
            pl.nth(0).alias("chrom"),
            pl.nth(1).cast(pl.Int64).alias("start"),
            pl.nth(2).cast(pl.Int64).alias("end"),
        )
        for p in paths
        if _has_rows(p)
    ]
    if not frames:
        return pl.DataFrame(
            schema={"chrom": pl.String, "start": pl.Int64, "end": pl.Int64}
        )
    return pl.concat(frames).filter(pl.col("chrom").is_in(list(sizes)))


def _has_rows(path):
    """Return whether a text file holds at least one line.

    Parameters
    ----------
    path : str
        File path.

    Returns
    -------
    bool
        True when the file is not empty.
    """
    with open(path) as f:
        return bool(f.readline())


def _chromosome(reader, chrom, length):
    """Read one chromosome of a bigWig.

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
        Signal per base as float32; 0 where the track has no data.
    """
    if chrom not in reader.chrom_sizes:
        return np.zeros(length, dtype=np.float32)
    return reader.read(chrom, [0], length)[0]


def _window_sums(signal, transcripts):
    """Sum the signal in [TSS-100, TSS+1000) of each transcript.

    Parameters
    ----------
    signal : numpy.ndarray
        Signal of the transcripts' chromosome and strand.
    transcripts : polars.DataFrame
        Transcripts on that chromosome and strand.

    Returns
    -------
    numpy.ndarray
        One sum per transcript.
    """
    cumulative = np.concatenate([[0.0], np.cumsum(signal, dtype=np.float64)])
    plus = transcripts["strand"] == "+"
    tss = np.where(plus, transcripts["start"], transcripts["end"] - 1)
    lo = np.where(plus, tss - TSS_UPSTREAM, tss - TSS_DOWNSTREAM + 1)
    hi = np.where(plus, tss + TSS_DOWNSTREAM, tss + TSS_UPSTREAM + 1)
    lo = np.clip(lo, 0, len(signal))
    hi = np.clip(hi, 0, len(signal))
    return cumulative[hi] - cumulative[lo]


def _region_sums(signal, transcripts):
    """Sum signal and length over the TSS regions and the gene body.

    Parameters
    ----------
    signal : numpy.ndarray
        Signal of the transcripts' chromosome and strand.
    transcripts : polars.DataFrame
        Chosen transcripts on that chromosome and strand.

    Returns
    -------
    numpy.ndarray
        TSS signal, TSS length, gene-body signal, gene-body length.
    """
    if transcripts.is_empty():
        return np.zeros(4)
    start = transcripts["start"].to_numpy()
    end = transcripts["end"].to_numpy()
    plus = transcripts["strand"].to_numpy() == "+"
    head = np.where(plus, start, end - END_REGION)
    tail = np.where(plus, end - END_REGION, start)
    tss = _mask_of(len(signal), head, head + END_REGION)
    body = (
        _mask_of(len(signal), start, end)
        & ~tss
        & ~_mask_of(len(signal), tail, tail + END_REGION)
    )
    return np.array(
        [
            signal[tss].sum(dtype=np.float64),
            tss.sum(),
            signal[body].sum(dtype=np.float64),
            body.sum(),
        ]
    )


def _mask(length, intervals):
    """Mark the bases that any interval covers.

    Parameters
    ----------
    length : int
        Chromosome length.
    intervals : polars.DataFrame
        Columns start and end, 0-based, half-open.

    Returns
    -------
    numpy.ndarray
        Boolean mask of the chromosome.
    """
    return _mask_of(length, intervals["start"].to_numpy(), intervals["end"].to_numpy())


def _mask_of(length, starts, ends):
    """Mark the bases in [starts, ends), clipped to the chromosome.

    Parameters
    ----------
    length : int
        Chromosome length.
    starts, ends : numpy.ndarray
        Interval bounds, 0-based, half-open.

    Returns
    -------
    numpy.ndarray
        Boolean mask of the chromosome.
    """
    edges = np.zeros(length + 1, dtype=np.int32)
    np.add.at(edges, np.clip(starts, 0, length), 1)
    np.add.at(edges, np.clip(ends, 0, length), -1)
    return np.cumsum(edges[:-1]) > 0


if __name__ == "__main__":
    main()
