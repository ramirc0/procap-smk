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
    """Read the tracks once per chromosome, then the windows QC needs.

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
    peaks = _peaks(args.peaks, sizes)

    totals = {"+": 0.0, "-": 0.0}
    in_peaks = 0.0
    for chrom, length in sizes.items():
        for strand, reader in readers.items():
            signal = _chromosome(reader, chrom, length)
            totals[strand] += float(signal.sum(dtype=np.float64))
            if peaks is not None:
                on = peaks.filter(pl.col("chrom") == chrom)
                in_peaks += sum(
                    float(signal[s:e].sum(dtype=np.float64))
                    for s, e in zip(on["start"], on["end"])
                )

    transcripts = _qc_transcripts(args.gtf, sizes)
    scores = pl.concat(
        on.with_columns(score=_tss_window_sums(readers[strand], on))
        for (strand,), on in transcripts.partition_by("strand", as_dict=True).items()
    )
    top = (
        scores.group_by("gene_id")
        .agg(pl.col("score").max())
        .sort("score", "gene_id", descending=[True, False])
        .head(max(1, int(TOP_FRACTION * scores["gene_id"].n_unique())))
    )
    # The chosen genes overlap no other gene, so summing per gene gives the
    # same totals as unions taken over all of them.
    sums = np.zeros(4)
    chosen = transcripts.join(top.select("gene_id"), on="gene_id")
    for gene in chosen.partition_by("gene_id"):
        sums += _gene_sums(readers[gene["strand"][0]], gene)
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
        Columns chrom, start, end of the merged peaks, so that overlapping
        calls count once.
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
    return (
        pl.concat(frames)
        .filter(pl.col("chrom").is_in(list(sizes)))
        .sort("chrom", "start")
        .with_columns(
            group=(pl.col("start") > pl.col("end").cum_max().shift(1))
            .fill_null(True)
            .cum_sum()
            .over("chrom")
        )
        .group_by("chrom", "group")
        .agg(pl.col("start").min(), pl.col("end").max())
        .drop("group")
    )


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


def _tss_window_sums(reader, transcripts):
    """Sum the signal in [TSS-100, TSS+1000) of each transcript, strand-aware.

    Parameters
    ----------
    reader : figwig.BigWigReader
        Track of the transcripts' strand.
    transcripts : polars.DataFrame
        Transcripts on that strand.

    Returns
    -------
    numpy.ndarray
        One sum per transcript.
    """
    plus = (transcripts["strand"] == "+").to_numpy()
    tss = np.where(plus, transcripts["start"], transcripts["end"] - 1)
    starts = np.where(plus, tss - TSS_UPSTREAM, tss - TSS_DOWNSTREAM + 1)
    windows = reader.read(
        transcripts["chrom"].to_numpy(),
        np.maximum(starts, 0),
        TSS_UPSTREAM + TSS_DOWNSTREAM,
    )
    return np.nansum(windows, axis=1, dtype=np.float64)


def _gene_sums(reader, gene):
    """Sum signal and length over one gene's TSS regions and gene body.

    Parameters
    ----------
    reader : figwig.BigWigReader
        Track of the gene's strand.
    gene : polars.DataFrame
        The gene's chosen transcripts.

    Returns
    -------
    numpy.ndarray
        TSS signal, TSS length, gene-body signal, gene-body length.
    """
    lo, hi = gene["start"].min(), gene["end"].max()
    signal = reader.read(gene["chrom"][0], [lo], hi - lo)[0]
    tss = np.zeros(hi - lo, dtype=bool)
    tts = np.zeros(hi - lo, dtype=bool)
    body = np.zeros(hi - lo, dtype=bool)
    for start, end, strand in gene.select("start", "end", "strand").iter_rows():
        start, end = start - lo, end - lo
        head, tail = (
            (start, end - END_REGION) if strand == "+" else (end - END_REGION, start)
        )
        body[start:end] = True
        tss[head : head + END_REGION] = True
        tts[tail : tail + END_REGION] = True
    body &= ~tss & ~tts
    return np.array(
        [
            signal[tss].sum(dtype=np.float64),
            tss.sum(),
            signal[body].sum(dtype=np.float64),
            body.sum(),
        ]
    )


if __name__ == "__main__":
    main()
