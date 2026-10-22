# Config, run sheet, unit maps, and shared helpers.

import re
from pathlib import Path

import polars as pl
from snakemake.exceptions import WorkflowError


configfile: "config/config.yaml"


RESULTS = config["outdir"]
RUN_ID = config.get("run_id", "default")
OUTDIR = f"{RESULTS}/{RUN_ID}"          # per-run
REFDIR = f"{RESULTS}/refs"              # run-independent
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"
REFERENCES = config["references"]
GENOMES = ["hg38_rDNA", "hg38_rDNA_dm6"]

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)


# The pixi plugin installs the env inside each job. Parallel SLURM jobs then race
# on a fresh or relocked env, so install it once before any job starts.
onstart:
    shell("pixi install --manifest-path {workflow.basedir:q}/envs/pixi.toml --all --locked")


_SHEET = Path(config["samples"])
_COLUMNS = ["run", "library", "experiment", "r1", "r2", "spikein"]
if not _SHEET.exists():
    raise WorkflowError(
        f"Run sheet not found: {_SHEET}\n"
        f"Columns: {', '.join(_COLUMNS)}; see config/samples.example.tsv."
    )

# All-string read so empty cells are "" (absent), not a type error.
_manifest = pl.read_csv(_SHEET, separator="\t", infer_schema_length=0, comment_prefix="#")
_missing = set(_COLUMNS) - set(_manifest.columns)
if _missing:
    raise WorkflowError(f"Run sheet {_SHEET} lacks columns: {sorted(_missing)}")


def _column(name):
    """Map each run to its value in sheet column `name`.

    Parameters
    ----------
    name : str
        Run sheet column.

    Returns
    -------
    dict of str to str or None
        Value per `run`, with empty cells as None.
    """
    return {row["run"]: (row[name] or None) for row in _manifest.iter_rows(named=True)}


RUNS = _manifest.get_column("run").to_list()
LIBRARY_OF = _column("library")
EXPERIMENT_OF = _column("experiment")
R1_OF = _column("r1")
R2_OF = _column("r2")
SPIKEIN_OF = _column("spikein")

RUNS_OF = {}
for _run in RUNS:
    RUNS_OF.setdefault(LIBRARY_OF[_run], []).append(_run)
LIBRARIES_OF = {}
for _library, _runs in RUNS_OF.items():
    LIBRARIES_OF.setdefault(EXPERIMENT_OF[_runs[0]], []).append(_library)
LIBRARIES = list(RUNS_OF)
EXPERIMENTS = list(LIBRARIES_OF)

if not RUNS:
    raise WorkflowError(f"Run sheet {_SHEET} has no rows.")
if len(set(RUNS)) != len(RUNS):
    raise WorkflowError(f"Run sheet {_SHEET} repeats a run.")
if set(RUNS) & set(LIBRARIES) or set(RUNS) & set(EXPERIMENTS) or set(LIBRARIES) & set(EXPERIMENTS):
    raise WorkflowError("Run, library and experiment IDs must not overlap; they share output directories.")
for _library, _runs in RUNS_OF.items():
    if len({EXPERIMENT_OF[r] for r in _runs}) > 1:
        raise WorkflowError(f"Library {_library} spans several experiments.")
    if len({SPIKEIN_OF[r] for r in _runs}) > 1:
        raise WorkflowError(f"Runs of library {_library} disagree on spikein; they are merged.")
if set(SPIKEIN_OF.values()) - {None, "dm6"}:
    raise WorkflowError("Column spikein takes dm6 or an empty cell.")

def _alternatives(names):
    """Return a regex that matches exactly one of `names`.

    Parameters
    ----------
    names : list of str
        Literal names.

    Returns
    -------
    str
        Alternation of the escaped names.
    """
    return "|".join(re.escape(n) for n in names)


wildcard_constraints:
    genome=_alternatives(GENOMES),
    run=_alternatives(RUNS),
    library=_alternatives(LIBRARIES),
    experiment=_alternatives(EXPERIMENTS),


def prefix(unit):
    """Return the output prefix of a run, library or experiment.

    Parameters
    ----------
    unit : str
        Run, library or experiment ID from the sheet.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<unit>/<unit>`.
    """
    return f"{OUTDIR}/{unit}/{unit}"


def genome_of(run):
    """Return the genome a run aligns to.

    Parameters
    ----------
    run : str
        Run ID from the sheet.

    Returns
    -------
    str
        `hg38_rDNA_dm6` for a dm6 spike-in run, else `hg38_rDNA`.
    """
    return "hg38_rDNA_dm6" if SPIKEIN_OF[run] == "dm6" else "hg38_rDNA"
