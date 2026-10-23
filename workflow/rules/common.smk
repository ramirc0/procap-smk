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
PEAK_TYPES = ["bidirectional", "divergent", "unidirectional"]

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


if not _manifest.height:
    raise WorkflowError(f"Run sheet {_SHEET} has no rows.")
if _manifest.select(pl.col("run", "library", "experiment").str.contains(r"\.").any()).row(0) != (False,) * 3:
    raise WorkflowError("Run, library and experiment IDs must not contain dots; dots join them in file names.")

# A library is `<experiment>.<library>`, a run `<experiment>.<library>.<run>`.
# These stems name every output file.
RUNS = []
R1_OF, R2_OF, SPIKEIN_OF, RUNS_OF = {}, {}, {}, {}
for _row in _manifest.iter_rows(named=True):
    _library = f"{_row['experiment']}.{_row['library']}"
    _run = f"{_library}.{_row['run']}"
    RUNS.append(_run)
    R1_OF[_run], R2_OF[_run] = _row["r1"], _row["r2"]
    SPIKEIN_OF[_run] = _row["spikein"] or None
    RUNS_OF.setdefault(_library, []).append(_run)
LIBRARIES = list(RUNS_OF)
LIBRARIES_OF = {}
for _library in LIBRARIES:
    LIBRARIES_OF.setdefault(_library.split(".")[0], []).append(_library)
EXPERIMENTS = list(LIBRARIES_OF)
REPLICATED = [e for e in EXPERIMENTS if len(LIBRARIES_OF[e]) > 1]

if len(set(RUNS)) != len(RUNS):
    raise WorkflowError(f"Run sheet {_SHEET} repeats a run.")
for _library, _runs in RUNS_OF.items():
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
    unit=_alternatives(LIBRARIES + EXPERIMENTS),


def out(rule, unit):
    """Return the output prefix of a rule for a run, library or experiment.

    Parameters
    ----------
    rule : str
        Rule that writes the output.
    unit : str
        Run or library stem, or experiment ID.

    Returns
    -------
    str
        Path `<outdir>/<run_id>/<experiment>/<rule>/<unit>`.
    """
    return f"{OUTDIR}/{unit.split('.')[0]}/{rule}/{unit}"


def genome_of(run):
    """Return the genome a run aligns to.

    Parameters
    ----------
    run : str
        Run stem.

    Returns
    -------
    str
        `hg38_rDNA_dm6` for a dm6 spike-in run, else `hg38_rDNA`.
    """
    return "hg38_rDNA_dm6" if SPIKEIN_OF[run] == "dm6" else "hg38_rDNA"


def libraries_in(unit):
    """Return the libraries a library or experiment track pools.

    Parameters
    ----------
    unit : str
        Library stem or experiment ID.

    Returns
    -------
    list of str
        The library itself, or every library of the experiment.
    """
    return [unit] if unit in RUNS_OF else LIBRARIES_OF[unit]


def peak_beds(experiment):
    """Return the PINTS peak BEDs of an experiment.

    Parameters
    ----------
    experiment : str
        Experiment ID from the sheet.

    Returns
    -------
    list of str
        Bidirectional, divergent and unidirectional peak BED paths.
    """
    return [f"{out('pints', experiment)}.{t}.bed" for t in PEAK_TYPES]
