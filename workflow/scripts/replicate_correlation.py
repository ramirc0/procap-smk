"""Pearson correlation between the libraries of one PRO-cap experiment.

As in Shah et al. (bioRxiv 10.1101/2025.09.24.676871), each track is summed
into 50 bp bins, plus and minus bins are pooled, and each library pair is
correlated on log10 counts over the bins non-zero in both.
"""

import argparse
import itertools

import numpy as np
import polars as pl
from figwig import BigWigReader


def build_parser():
    """Build the replicate_correlation.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the replicate_correlation.py flags.
    """
    p = argparse.ArgumentParser(description="Correlate the libraries of an experiment.")
    p.add_argument("--experiment", required=True, help="Experiment ID.")
    p.add_argument("--libraries", nargs="+", required=True, help="Library IDs.")
    p.add_argument(
        "--plus", nargs="+", required=True, help="Plus-strand bigWigs, one per library."
    )
    p.add_argument(
        "--minus",
        nargs="+",
        required=True,
        help="Minus-strand bigWigs, one per library.",
    )
    p.add_argument("--chrom_sizes", required=True, help="Chromosomes to bin.")
    p.add_argument("--bin_size", type=int, default=50, help="Bin width in bp.")
    p.add_argument("-o", "--output", required=True, help="TSV with one row per pair.")
    return p


def main(argv=None):
    """Bin every library, then correlate each pair.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    if not len(args.libraries) == len(args.plus) == len(args.minus):
        raise SystemExit("Give one plus and one minus bigWig per library.")
    sizes = pl.read_csv(args.chrom_sizes, separator="\t", has_header=False).rows()
    bins = {
        lib: np.concatenate(
            [
                _binned(BigWigReader(path), sizes, args.bin_size)
                for path in (plus, minus)
            ]
        )
        for lib, plus, minus in zip(args.libraries, args.plus, args.minus)
    }
    rows = []
    for a, b in itertools.combinations(args.libraries, 2):
        both = (bins[a] > 0) & (bins[b] > 0)
        r = (
            np.corrcoef(np.log10(bins[a][both]), np.log10(bins[b][both]))[0, 1]
            if both.sum() > 1
            else np.nan
        )
        rows.append((args.experiment, a, b, int(both.sum()), float(r)))
    pl.DataFrame(
        rows,
        schema=["experiment", "library_a", "library_b", "shared_bins", "pearson_r"],
        orient="row",
    ).write_csv(args.output, separator="\t")


def _binned(reader, sizes, bin_size):
    """Sum a track into fixed bins over every chromosome.

    Parameters
    ----------
    reader : figwig.BigWigReader
        Open track.
    sizes : list of tuple of (str, int)
        Chromosomes and lengths, in order.
    bin_size : int
        Bin width in bp.

    Returns
    -------
    numpy.ndarray
        Bin sums, chromosome after chromosome.
    """
    out = []
    for chrom, length in sizes:
        n = -(-length // bin_size)
        if chrom not in reader.chrom_sizes:
            out.append(np.zeros(n))
            continue
        signal = np.zeros(n * bin_size, dtype=np.float32)
        signal[:length] = reader.read(chrom, [0], length)[0]
        out.append(signal.reshape(n, bin_size).sum(axis=1, dtype=np.float64))
    return np.concatenate(out)


if __name__ == "__main__":
    main()
