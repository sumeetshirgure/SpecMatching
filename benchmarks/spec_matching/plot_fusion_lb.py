#!/usr/bin/env python3
# Copyright 2026 SpecMatching contributors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Draws the divide-and-conquer latencies `fusion_lb_profiler` logs.

    python benchmarks/spec_matching/plot_fusion_lb.py path/to/run_dir [more_run_dirs ...]

## What a run directory is, and why this reader takes many files

`load_balancer_profiler` wrote one `shots.csv` for a whole campaign. `fusion_lb_profiler` writes
**one CSV per cell** — `shots_d{d}_p{p}_T{T}_k{k}.csv` — beside a shared `agg.json` and `run.log`.
So a run directory is a *set* of files, and this script globs them, reads the cell identity off each
file's leading `# d=..,p=..,T=..,k=..` comment line, and sorts the set into `(d, p, T, k)` order
before reading a byte of data. Several run directories can be given at once and are processed one
after another, each producing its own figures and its own table.

The identity is read from the comment line rather than from the filename because the filename is a
`%.10g` rendering of the same numbers and a reader that parses it is a reader that can be confused by
`p=1e-05`. The filename is the fallback, and a file that has neither is skipped by name.

Rows are streamed in blocks, so a million-shot cell costs the block and not the campaign; what is
resident is the `|k|` cells of one `(d, p, T)` group, which is what the small multiples need.

## The two series, and why only two

Per shot the profiler logs three latencies, and `run.log` states the definitions:

    sparse_k = the makespan of the simulated k-core schedule (PRE, dispatches, leaves, fusions,
               extractions, COMBINE)
    fallback = stock exact decode on G, from the full syndrome, on its own core
    system   = escalated ? fallback : min(fallback, sparse_k)

The histograms draw **`fallback` and `sparse_k`**. `system` is a per-shot `min` of the two, so its
curve would sit on top of `sparse_k` everywhere but the escalating tail; it is a dashed reference
rule and it is in every table.

Two hues, and that is a checked decision rather than a stylistic one. The corpus palette's slot 0
(blue) against slot 1 (orange) is ΔE 24.7 under protanopia in the light theme and 26.8 in the dark
one, against a floor of 8; slot 1 against slot 2 (green) is 4.4, which no amount of direct labelling
rescues. So every figure here carries at most two hues, and anything needing more categories is a
small multiple. That is why the `k` sweep is a row of panels and not `|k|` curves on one axis, and
why every breakdown figure is single-hue with row labels carrying the identity.

## Escalating shots are not in the `sparse_k` series

This is the one place this script must not copy its predecessor. Under the profiler's deviation 5 an
escalating shot **stops at the job that truncated** and its `ticks_sparse_k` is written as `0` — the
schedule never finished, so there is no makespan to report. Feeding those zeroes to a histogram would
put a spike at the origin and drag the mean down, and the drag would grow with the escalation rate,
which is exactly the quantity a `T` sweep moves.

So `sparse_k` is drawn and averaged over **non-escalating shots only**, `fallback` and `system` over
all of them, and every caption and table says which is which. `agg.json` does the same thing for the
same reason, so the two agree.

## The figures

  * **one histogram per `(d, p, T, k)` cell** — `fallback` against `sparse_k`, means direct-labelled,
    the `system` mean as a reference rule;
  * **one `k`-sweep per `(d, p, T)`** — the same two histograms as small multiples, one panel per
    `k`, on shared bins and shared axes. `fallback` is the internal control: it is the same shots
    decoded the same way in every panel, so it should *not* move with `k`;
  * **one critical-path figure per `(d, p, T)`** — two rows of panels, one column per `k`. The top
    row is where the makespan goes (manager, build, leaf, fuse, extract, idle); the bottom row is
    what the manager's own share is made of (the three preprocessing regions, the dispatches, the
    combine). Linear axes, not log: a bar encodes length from zero, and a log axis would inflate
    `combine` into looking comparable with the work;
  * **one work figure per `(d, p, T)`** — solver time summed over the `k` cores against solver time
    on the critical path, as `k` grows. This is the question the binary exists to answer: cutting an
    oversized component is worth doing exactly insofar as it moves the second bar down while the
    first stays put;
  * **one structure figure per `(d, p, T)`** — the cut itself, from `agg.json`'s histograms: how many
    pieces, how deep the fusion tree, how big the largest fusion boundary, how many regions a fusion
    releases. Rows are the statistics, columns are `k`.

## Two identities, checked rather than assumed

`run.log` states both, and a figure drawn from data that no longer satisfies them is a figure
captioned with the wrong model. The table checks them and says so when they fail:

    cp_manager + cp_build + cp_leaf + cp_fuse + cp_extract + cp_idle == sparse_k mean
    uf + tree + scatter                                             == pre mean

`cp_idle` is additionally expected to be **exactly zero**: the profiler's critical path is walked
back through whatever determined each start, so every link is tight by construction. A non-zero
`cp_idle` means the schedule and its own critical path disagree, which is a bug rather than a
measurement, and the table flags it.

## Overhead, and why the correction is per-shot here

Rows are raw, uncorrected tick deltas; `run.log` carries `timer_overhead_ns` and says the analysis is
where the subtraction happens. This is that analysis.

A bracketed region costs about one clock read, so a column carries one `timer_overhead_ns` per
interval summed into it. `fallback` is one interval. `sparse_k` is not a fixed count — it is a
makespan over however many jobs happened to land on the critical path — but the profiler logs that
count per shot as `cp_len`, plus the two interior reads inside `PRE` that split it into `uf`, `tree`
and `scatter`. So `--subtract-overhead` takes `(cp_len + 2)` reads off `sparse_k` **per shot** and
one off `fallback`, which is a sharper correction than any flat one and needs no assumption about the
shape of the schedule.

That is a **lower bound** by one read per leaf on the path, and deliberately so. A `LEAF` is now two
timed regions — the graph build and the solve — and so carries one interior read of its own, but how
many leaves a given shot's path ran through is not a logged quantity; `cp_len` counts jobs, not
regions. Usually it is one leaf, and on the `rdpmc` backend one read is tens of nanoseconds against a
makespan in the tens of microseconds. Inventing a count would be worse than under-correcting by a
known amount and saying so, which is the same rule the `cp_*` means are left raw under.

The means in `agg.json` are corrected at one interval per module. The `cp_*` means are left **raw**,
because the split of `cp_len` across the six buckets is not logged and inventing one would be worse
than saying so; the table prints the mean `cp_len` beside them so a reader can bound it. The two
identities above are therefore always checked on raw numbers.

## Which order the edges arrived in

The profiler's `--sort` — the hardware sort that hands the union-find `H`'s edges in non-decreasing
weight — is **off by default**, so "a cut falls on a component's heaviest links" is a property of the
run and not of the binary. It changes which edges are refused, and with them the piece count, the
fusion tree, the largest fusion boundary, how many regions a fusion releases and the makespan, which
is every panel of the structure figure and most of the critical-path one. So it is read off
`run.log`'s `edge_order=` line (`agg.json`'s `sort_edges` is the fallback) and printed in the table
header and on both of those figures; no caption in this script asserts it. It also decides whether the
`cut` column of the gate table means anything: the two cut-on-heaviest checks are claims about a
sorted input and the profiler skips them without `--sort`, so the table says "not checked" rather than
letting a `0` read as a pass.

## What travels with every number

`run.log` states the model and a caveat, and they are not optional context: one shared solver
instance holds the whole shot, **no cost is modelled for moving data between cores**, every unit of
work runs serially on one thread and is placed on a simulated schedule, and a unit of work executed
straight after another sees a warmer cache than a real core would. The captions repeat the caveat and
the table repeats the model. A figure from this script without them is being read as something it is
not.
"""

import argparse
import glob
import itertools
import json
import math
import os
import re
import sys
import textwrap

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

# The corpus palette, shared with the other plotters here so a reader moving between their output
# sees one visual language.
#
# Only slots 0 and 1 are used for series. Checked with the dataviz validator against each theme's own
# surface rather than trusted: light passes at ΔE 24.7 protan / 33.6 normal, dark at 26.8 / 31.8, and
# both clear the lightness band, the chroma floor and 3:1 contrast. Slot 1 against slot 2 (green) is
# ΔE 4.4 protan, which is why there is no third series anywhere in this file.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#dcdcd8",
        # (fallback on G, the fused k-core path).
        "series": ("#2a78d6", "#eb6834"),
    },
    "dark": {
        "surface": "#1a1a19",
        "text_primary": "#ffffff",
        "text_secondary": "#c3c2b7",
        "grid": "#3a3a37",
        "series": ("#3987e5", "#d95926"),
    },
}

FALLBACK_LABEL = "fallback: stock decode on G"
SPARSE_LABEL = "sparse_k: fused k-core makespan on H"

# Rows are parsed in blocks of this many, so the transient cost of a read is set by the block and not
# by the length of the campaign.
CHUNK_ROWS = 100_000

# The only columns any series needs. Everything else the figures show is a mean, and `agg.json`
# already has those accumulated over the whole cell in 128-bit integers, which is a better number
# than anything re-summed here in floats.
#
# `cp_len` is read for the per-shot overhead correction and for the table; `ticks_system` is read
# even though it is recomputed here, so the recomputation can be checked rather than trusted.
NEEDED = (
    "escalated",
    "cp_len",
    "ticks_sparse_k",
    "ticks_fallback",
    "ticks_system",
)

# Timed intervals summed into each column. `fallback` is one bracketed region. `sparse_k` is a
# makespan over the jobs on the critical path — `cp_len` of them — plus the two interior clock reads
# inside `PRE` that split it into its three sub-regions. See the module docstring.
INTERVALS_FALLBACK = 1
PRE_INTERIOR_READS = 2
INTERVALS_MODULE = 1

# Above this share of clamped columns, `--subtract-overhead` is subtracting more than was measured
# and every corrected number is an artefact of the clamp. Said loudly rather than left in a footer.
CLAMP_ALARM_PERCENT = 1.0

# The makespan's own decomposition, in the order the schedule composes it. These six sum to the
# `sparse_k` mean.
#
# `build` is the leaf's graph build and `leaf` is its solve: since the graph build moved off the
# manager and into the leaves, a `LEAF` job is two consecutive timed regions on one core and the
# profiler charges the critical path for each separately.
CRITICAL_PATH = ("manager", "build", "leaf", "fuse", "extract", "idle")

# The manager core's work. The first three sum to `pre`; `dispatch_total` and `combine` are the
# manager's other two jobs and are shown beside them rather than in a second figure.
#
# There is no `bucket`, in either edge order: the design's bucketing pass is gone, and where a sort
# does run (`--sort`) it is the hardware's, stood in for outside every timed region. `run.log` records
# which order the run used as `edge_order=`, and `edge_order_note` is what puts it on the figure —
# nothing here may assert it, because `--sort` is off by default. A run old enough to have the column
# reads as a missing bar here, which is the right way round — a reader that silently averaged it into
# `pre` would be wrong.
MANAGER_MODULES = ("uf", "tree", "scatter", "dispatch_total", "combine")
PRE_MODULES = ("uf", "tree", "scatter")

# Solver work summed over the `k` cores, against the part of it that lands on the critical path.
# The build is solver-core work like the solve and the fusion, so it is in both sides of the
# comparison; leaving it out of the left bar would make the cut look free.
WORK_TOTAL = ("build_sum", "solve_sum", "fuse_sum")
WORK_PATH = ("build", "leaf", "fuse")

# The cut's own structure, as `agg.json` histograms them.
STRUCTURE = (
    ("n_pieces_hist", "pieces per shot"),
    ("tree_depth_hist", "fusion tree depth"),
    ("max_fusion_boundary_hist", "largest fusion boundary"),
    ("n_released_hist", "regions released"),
)

CELL_FILE_RE = re.compile(r"^shots_d.*_p.*_T.*_k.*\.csv$")
# `# d=7,p=0.003,T=1.5,k=8,alpha=1,seed=...` — the profiler's own first line, and the authority on
# what cell a file holds.
COMMENT_KEY_RE = re.compile(r"([A-Za-z_]+)=([^,\s]+)")


class Run:
    """One profiler output directory: its `run.log` metadata, its `agg.json` cells, and the list of
    per-cell CSV files it holds.

    Deliberately not the rows. Those stay on disk and are streamed by `stream_groups`; a `Run` is
    small enough that `main` can hold one per directory for the whole invocation.
    """

    def __init__(self, path, meta, agg, cell_files):
        self.path = path
        self.meta = meta
        self.agg = agg
        # `(d, p, T, k) -> path`, sorted so that a `(d, p, T)` group's files are contiguous and `k`
        # ascends within it.
        self.cell_files = cell_files
        # Set when the caller named one cell file rather than the directory; see `load`.
        self.only = None
        self.cells = {
            (float(cell["d"]), float(cell["p"]), float(cell["T"]), float(cell["k"])): cell
            for cell in agg.get("cells", ())
        }

    @property
    def name(self):
        return os.path.basename(os.path.normpath(self.path)) or self.path

    @property
    def timer(self):
        return self.meta.get("hires_timer_name", "?")

    @property
    def thread_scoped(self):
        return self.meta.get("timer_is_thread_scoped") == "1"

    @property
    def overhead_ns(self):
        """`timer_overhead_ns=9.2114 (10^6 back-to-back ...)` — the number, not the parenthetical."""
        text = self.meta.get("timer_overhead_ns", "")
        try:
            return float(text.split()[0])
        except (IndexError, ValueError):
            return float("nan")

    @property
    def alpha(self):
        try:
            return float(self.meta.get("alpha", "nan").split()[0])
        except (IndexError, ValueError):
            return float("nan")

    @property
    def verified(self):
        return self.meta.get("verify", "0").startswith("1")

    @property
    def edge_order(self):
        """`sorted_by_hardware` or `input_order`, off `run.log`'s `edge_order=` line.

        The profiler's `--sort` is off by default, so whether the union-find saw `H`'s edges in
        non-decreasing weight is a property of the *run* and no longer one of the binary. It is the
        difference between a cut on a component's heaviest links and a cut wherever the `ceil(n/k)`
        threshold happened to bite, and it moves the structure figures and the makespan with it, so
        every caption that used to assert the sorted order reads this instead.

        `agg.json`'s `sort_edges` is the fallback for a log whose line was lost. `None` means neither
        said, which is what a run older than the flag looks like if its log is gone — and is reported
        as unknown rather than guessed at.
        """
        stated = self.meta.get("edge_order", "").split()
        if stated:
            return stated[0]
        flag = self.agg.get("sort_edges")
        if flag is True:
            return "sorted_by_hardware"
        if flag is False:
            return "input_order"
        return None

    @property
    def sorted_edges(self):
        """Whether this run sorted its edge list by weight; `None` when the run does not say."""
        order = self.edge_order
        return None if order is None else order.startswith("sorted")


def read_run_log(path):
    """The top-level `key=value` lines of `run.log`.

    The indented `corpus` / `cell` lines below them are per-cell provenance the figures do not read —
    `agg.json` carries the same cells in a form that does not need parsing — so they are skipped by
    the leading-whitespace test rather than half-parsed. Values may themselves contain `=`, which is
    why this partitions on the first one instead of splitting.
    """
    meta = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line[:1].isspace():
                continue
            key, sep, value = line.partition("=")
            if sep:
                meta[key.strip()] = value.strip()
    return meta


def cell_key_of(path):
    """`(d, p, T, k)` for one per-cell CSV, from its leading comment line.

    The comment line is the authority: it is written by the profiler from the same values that named
    the file, and unlike the name it is unambiguous about `p=1e-05`. The name is the fallback for a
    file whose first line was lost, and a file with neither is not a cell file.
    """
    try:
        with open(path) as handle:
            first = handle.readline()
    except OSError:
        return None
    if first.startswith("#"):
        fields = dict(COMMENT_KEY_RE.findall(first[1:]))
        try:
            return (float(fields["d"]), float(fields["p"]), float(fields["T"]), float(fields["k"]))
        except (KeyError, ValueError):
            pass
    match = re.match(r"^shots_d([^_]+)_p([^_]+)_T([^_]+)_k([^.]+)\.csv$", os.path.basename(path))
    if not match:
        return None
    try:
        return tuple(float(value) for value in match.groups())
    except ValueError:
        return None


def load(path):
    """Reads one run directory's metadata and its list of cell files. Rows are left on disk.

    Everything that would make a figure wrong is settled here, before a figure is started: that the
    directory holds the three kinds of file, that at least one cell file is readable and names its
    cell, that its header has the columns the series are built from, and that `agg.json` parses. A
    malformed run is reported by name like any other bad input.
    """
    # A single cell file means that cell, not its whole directory: the `agg.json` and `run.log`
    # beside it are still needed (they carry the means and the model), but a reader who names one
    # file and is handed the other twenty-three has been surprised by their own tooling.
    only = None
    if os.path.isfile(path) and CELL_FILE_RE.match(os.path.basename(path)):
        only = cell_key_of(path)
        path = os.path.dirname(path) or "."
    if not os.path.isdir(path):
        raise ValueError("not a directory (expected a fusion_lb_profiler --out dir)")
    log_path = os.path.join(path, "run.log")
    agg_path = os.path.join(path, "agg.json")
    for needed in (log_path, agg_path):
        if not os.path.isfile(needed):
            raise ValueError(f"{path}: no {os.path.basename(needed)}")

    found = []
    for candidate in sorted(glob.glob(os.path.join(path, "shots_d*_p*_T*_k*.csv"))):
        key = cell_key_of(candidate)
        if key is None:
            print(f"  {candidate}: no cell identity in the comment line or the name; skipped", file=sys.stderr)
            continue
        found.append((key, candidate))
    if not found:
        raise ValueError(f"{path}: no shots_d*_p*_T*_k*.csv cell files")
    # `(d, p, T, k)` order with `k` innermost, so a group's files are contiguous and ascend in `k`.
    found.sort(key=lambda item: item[0])

    columns = header_of(found[0][1])
    missing = [name for name in NEEDED if name not in columns]
    if missing:
        raise ValueError(f"{found[0][1]}: header is missing {', '.join(missing)}")

    meta = read_run_log(log_path)
    with open(agg_path) as handle:
        agg = json.load(handle)
    if agg.get("binary") not in (None, "fusion_lb_profiler"):
        raise ValueError(f"{agg_path}: binary is {agg['binary']}; this script reads fusion_lb_profiler")
    run = Run(path, meta, agg, found)
    run.columns = columns
    run.only = only
    return run


def header_of(path):
    """The column names of one cell file, skipping its leading `#` comment line."""
    with open(path) as handle:
        line = handle.readline()
        if line.startswith("#"):
            line = handle.readline()
        if not line.strip():
            raise ValueError(f"{path}: no column header")
        return [name.strip() for name in line.strip().split(",")]


def parse_number_list(text, cast):
    """`--k 3,5` / `--T 1.5,2` as a set, or None when the flag was not given."""
    if text is None:
        return None
    values = set()
    for item in text.split(","):
        item = item.strip()
        if item:
            values.add(cast(item))
    return values or None


class Cell:
    """One `(d, p, T, k)` cell's drawable series, accumulated a block at a time.

    Holds three `float64` arrays in microseconds and a few scalars. It never holds the module columns
    or the structure counters: those are means and distributions, and `agg.json` already has them
    over the whole cell in 128-bit integers.

    `sparse` is shorter than `fallback` whenever the cell escalated anything — see the module
    docstring. That asymmetry is the point, not an accident, and `escalated` is what explains it.
    """

    def __init__(self, key, agg):
        self.key = key
        self.d, self.p, self.T, self.k = key
        self.agg = agg
        self._fallback = []
        self._sparse = []
        self._system = []
        self.escalated = 0
        self.cp_len_sum = 0.0
        self.shots = 0
        # Rows where the recomputed `system` disagreed with the column the profiler logged. Only
        # meaningful on the raw path; the corrected path recomputes from corrected inputs by design.
        self.system_mismatch = 0
        # Corrected values that came out below zero and were clamped. A negative latency is not a
        # thing, but silently clamping one is worse than saying how many.
        self.clamped = 0

    def add(self, fallback, sparse, system, escalated, cp_len):
        self._fallback.append(fallback)
        self._sparse.append(sparse)
        self._system.append(system)
        self.escalated += int(np.count_nonzero(escalated))
        self.cp_len_sum += float(cp_len.sum())
        self.shots += int(fallback.size)

    def finish(self):
        empty = np.empty(0, dtype=np.float64)
        self.fallback = np.concatenate(self._fallback) if self._fallback else empty
        self.sparse = np.concatenate(self._sparse) if self._sparse else empty
        self.system = np.concatenate(self._system) if self._system else empty
        self._fallback.clear()
        self._sparse.clear()
        self._system.clear()
        return self

    def release(self):
        self.fallback = self.sparse = self.system = np.empty(0)

    @property
    def escalation_rate(self):
        return self.escalated / self.shots if self.shots else float("nan")

    @property
    def cp_len(self):
        return self.cp_len_sum / self.shots if self.shots else float("nan")

    @property
    def title(self):
        return f"d = {self.d:g},  p = {self.p:g},  T = {self.T:g},  k = {self.k:g}"

    @property
    def slug(self):
        return f"d{self.d:g}_p{self.p:g}_T{self.T:g}_k{self.k:g}"

    def module(self, name, overhead=0.0):
        """One `agg.json` module mean, in raw nanoseconds, optionally overhead-corrected."""
        means = (self.agg or {}).get("module_means") or {}
        if name not in means:
            return float("nan")
        return max(float(means[name]) - INTERVALS_MODULE * overhead, 0.0)

    def path_mean(self, name):
        """One `agg.json` critical-path mean, in raw nanoseconds. Never overhead-corrected; see the
        module docstring."""
        means = (self.agg or {}).get("critical_path_means") or {}
        return float(means.get(name, float("nan")))

    def histogram(self, name):
        """One `agg.json` structure histogram as `(values, counts)`, ascending."""
        raw = (self.agg or {}).get(name) or {}
        if not raw:
            return np.empty(0), np.empty(0)
        items = sorted((int(key), int(value)) for key, value in raw.items())
        return (
            np.array([item[0] for item in items], dtype=np.float64),
            np.array([item[1] for item in items], dtype=np.float64),
        )

    def histogram_mean(self, name):
        values, counts = self.histogram(name)
        total = counts.sum()
        return float((values * counts).sum() / total) if total else float("nan")


def group_slug(key):
    d, p, T = key
    return f"d{d:g}_p{p:g}_T{T:g}"


def group_title(key):
    d, p, T = key
    return f"d = {d:g},  p = {p:g},  T = {T:g}"


def read_cell(run, key, path, args):
    """One cell file, streamed in blocks into a `Cell`.

    The correction, when asked for, is applied here — while the block is still a block — so a
    corrected campaign costs one extra vector op per block rather than a second pass.
    """
    at = {name: position for position, name in enumerate(NEEDED)}
    columns = header_of(path)
    usecols = tuple(columns.index(name) for name in NEEDED)
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    cell = Cell(key, run.cells.get(key))

    with open(path) as handle:
        line = handle.readline()
        if line.startswith("#"):
            handle.readline()
        while True:
            text = list(itertools.islice(handle, CHUNK_ROWS))
            if not text:
                break
            block = np.loadtxt(text, delimiter=",", usecols=usecols, dtype=np.float64, ndmin=2)
            del text
            if not block.size:
                continue
            escalated = block[:, at["escalated"]] != 0
            cp_len = block[:, at["cp_len"]]
            fallback = block[:, at["ticks_fallback"]]
            sparse = block[:, at["ticks_sparse_k"]]
            if overhead:
                # One read per bracketed region: one for `fallback`, and for `sparse_k` however many
                # jobs landed on this shot's critical path plus the three interior reads inside
                # `PRE`. See the module docstring.
                fallback = fallback - INTERVALS_FALLBACK * overhead
                sparse = sparse - (cp_len + PRE_INTERIOR_READS) * overhead
                below = int(np.count_nonzero(fallback < 0) + np.count_nonzero(sparse[~escalated] < 0))
                if below:
                    cell.clamped += below
                    np.clip(fallback, 0.0, None, out=fallback)
                    np.clip(sparse, 0.0, None, out=sparse)
            # The definition `run.log` states, recomputed here rather than read, so that the
            # corrected path and the raw path produce it the same way.
            system = np.where(escalated, fallback, np.minimum(fallback, sparse))
            if not overhead:
                cell.system_mismatch += int(np.count_nonzero(system != block[:, at["ticks_system"]]))
            # An escalating shot has no makespan — it stopped at the job that truncated, and the
            # profiler writes `ticks_sparse_k = 0`. Those rows leave the sparse series here rather
            # than being drawn as a spike at the origin.
            cell.add(
                fallback / 1000.0,
                sparse[~escalated] / 1000.0,
                system / 1000.0,
                escalated,
                cp_len,
            )
            del block
    return cell.finish()


def stream_groups(run, args):
    """Yields `((d, p, T), [Cell, ...])`, one group at a time, `k` ascending within it.

    `run.cell_files` is already sorted by `(d, p, T, k)`, so a group is a run of adjacent entries and
    can be completed, handed to the figures and released before the next one is opened. What is
    resident is `|k|` cells rather than the campaign.
    """
    current_key = None
    pending = []
    for key, path in run.cell_files:
        if run.only is not None and key != run.only:
            continue
        if args.d is not None and key[0] not in args.d:
            continue
        if args.p is not None and key[1] not in args.p:
            continue
        if args.T is not None and key[2] not in args.T:
            continue
        if args.k is not None and key[3] not in args.k:
            continue
        group = key[:3]
        if group != current_key:
            if pending:
                yield current_key, pending
            current_key = group
            pending = []
        cell = read_cell(run, key, path, args)
        if cell.shots:
            pending.append(cell)
    if pending:
        yield current_key, pending


def wrap_caption(text, width=135):
    """Re-wraps a caption to a character width, keeping the author's own line breaks.

    A caption that runs off the right edge is a caption that was not read, and these say things — the
    model, the warm-cache caveat, which shots are not in which series, which order the edges arrived
    in — that the figure is not honest without.

    The width is the caller's, and `caption_block` derives it from the figure it is about to write
    on: these range from four inches across (a one-cell run's small multiples) to fourteen (a four-k
    sweep), and one constant cannot be right for both. 135 is the old wide default, kept for a direct
    caller that has no figure in hand.
    """
    return "\n".join(
        textwrap.fill(line, width=width, break_long_words=False) if line else ""
        for line in text.split("\n")
    )


def percentile(values, fraction):
    """Only ever used to choose where the x axis stops; nothing reported here is a percentile.

    A partial sort, not a full one: the axis needs one order statistic. The percentiles that *are*
    reported come from `agg.json`, where the profiler computed them over the raw ticks.
    """
    if values.size == 0:
        return float("nan")
    index = min(int(fraction * (values.size - 1)), values.size - 1)
    return float(np.partition(values, index)[index])


def mean_of(values):
    return float(values.mean()) if values.size else float("nan")


def stats_of(values):
    return {
        "n": int(values.size),
        "mean": mean_of(values),
        "max": float(values.max()) if values.size else float("nan"),
    }


def bin_edges(sides, x_max, bins, log_x):
    """Shared edges for every series and every panel — histograms on different bins are not
    comparable, and the `k` sweep exists precisely to be compared across panels."""
    if log_x:
        lows = [float(side[side > 0].min()) for side in sides if np.any(side > 0)]
    else:
        lows = [float(side.min()) for side in sides if side.size]
    low = max(min(lows), 1e-3) if lows else 1e-3
    high = max(x_max, low * (1.0 + 1e-6))
    if log_x:
        return np.logspace(math.log10(low), math.log10(high), bins + 1)
    return np.linspace(low, high, bins + 1)


def axis_limits(cells, args):
    """The x range and bins shared by every panel of a group, from the two drawn series only.

    `system` does not widen the axis: it is bounded above by `fallback` on every shot, so an axis
    that covers the two inputs covers it too.
    """
    sides = []
    for cell in cells:
        sides.extend((cell.fallback, cell.sparse))
    pooled = np.concatenate([side for side in sides if side.size]) if sides else np.empty(0)
    x_max = max(
        percentile(pooled, args.x_max_percentile / 100.0),
        1.15 * max((mean_of(side) for side in sides if side.size), default=1.0),
    )
    del pooled
    return x_max, bin_edges(tuple(sides), x_max, args.bins, args.log_x)


def style_axes(ax, theme, log_x=False, log_y=False, grid_axis="y"):
    """The recessive frame every panel wears: grid behind the data, two spines, muted ticks."""
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")
    ax.grid(axis=grid_axis, color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)


def draw_histograms(
    ax,
    cell,
    edges,
    theme,
    log_x=False,
    label=True,
    series_labels=(FALLBACK_LABEL, SPARSE_LABEL),
    show_system=True,
    fontsize=9,
    linewidth=2.0,
    unit="us",
):
    """The two series into one axes, and the `system` mean as a reference rule.

    Binned once and drawn twice per series. `ax.hist` would re-bin the whole array for each pass and
    keep a copy inside the axes; past a few hundred thousand shots that is the difference between a
    figure and a swap storm, and the counts are identical either way.

    Filled at low alpha with a 2px edge, so the overlap is legible either way round.

    `series_labels`, `show_system`, `fontsize` and `linewidth` exist for `--paper`, which draws the
    same figure at column width with short labels and leaves the explanation to the paper's caption.
    """
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    for values, colour, series_label in (
        (cell.fallback, theme["series"][0], series_labels[0]),
        (cell.sparse, theme["series"][1], series_labels[1]),
    ):
        counts, _ = np.histogram(values, bins=edges)
        ax.stairs(counts, edges, color=colour, alpha=0.45, fill=True, label=series_label, zorder=2)
        ax.stairs(counts, edges, color=colour, linewidth=linewidth, zorder=3)

    def side_of(value):
        """Which side of its rule a label goes, as `(ha, x offset)`.

        A label sits to the right of its line by default and flips left once the line is past ~60% of
        the axis, where "to the right" runs it into the legend box or off the edge. The means sit
        close together at most grid points, so this is decided per label rather than once.
        """
        low, high = (math.log10(edges[0]), math.log10(edges[-1])) if log_x else (edges[0], edges[-1])
        at = (math.log10(value) if log_x else value) if value > 0 else low
        fraction = (at - low) / (high - low) if high > low else 0.0
        return ("right", -6) if fraction > 0.6 else ("left", 5)

    # Selective direct labels: the two means the figure is about, and the derived `system` mean as a
    # rule in text ink rather than a third series colour — it is a reference, not an identity.
    annotations = (
        (mean_of(cell.fallback), theme["series"][0], theme["text_primary"], 0.97, "mean {:,.2f} " + unit),
        (mean_of(cell.sparse), theme["series"][1], theme["text_primary"], 0.88, "mean {:,.2f} " + unit),
        (
            mean_of(cell.system),
            theme["text_secondary"],
            theme["text_secondary"],
            0.79,
            "system mean {:,.2f} " + unit,
        ),
    )
    for index, (mean, rule_colour, text_colour, height, template) in enumerate(annotations):
        if not math.isfinite(mean):
            continue
        system = index == 2
        if system and not show_system:
            continue
        ax.axvline(
            mean,
            color=rule_colour,
            linewidth=1.2 if system else 1.5,
            linestyle=(0, (1, 2)) if system else (0, (4, 3)),
            zorder=4,
        )
        if not label:
            continue
        ha, offset = side_of(mean)
        ax.annotate(
            template.format(mean),
            xy=(mean, height),
            xycoords=("data", "axes fraction"),
            ha=ha,
            va="top",
            fontsize=fontsize,
            color=text_colour,
            xytext=(offset, 0),
            textcoords="offset points",
            zorder=5,
        )


def model_caveat(run):
    """The model and the caveat that travel with every number, in one short line.

    Taken from `run.log` rather than restated here, so a run written under a changed model cannot be
    captioned with this script's memory of the old one.
    """
    parts = []
    if run.meta.get("fusion_k_core_model", "").startswith("1"):
        parts.append(
            "fusion_k_core_model: one shared instance holds the whole shot; every unit of work runs"
            " serially on one thread and is placed on a simulated k-core schedule"
        )
    if run.meta.get("shared_instance", "").startswith("1"):
        parts.append("no cost is modelled for moving data between cores")
    if run.meta.get("warm_cache_caveat", "").startswith("1"):
        parts.append("a unit of work run straight after another sees a warmer cache than a real core would")
    return "  ·  ".join(parts)


def edge_order_note(run):
    """Which order the union-find visited the edges in, in one line.

    Taken from `run.log` rather than restated here, for the same reason as `model_caveat`: `--sort`
    is off by default, so a caption that asserts the sorted order is a caption that is wrong about
    most runs. The cut, the fusion tree and every structure panel depend on this, which is why it
    goes on the figures and not only in the header.
    """
    sorted_edges = run.sorted_edges
    if sorted_edges is None:
        return "edge order: NOT STATED by this run — which links the cut fell on is unknown"
    if sorted_edges:
        return (
            "edge order: sorted by weight (--sort), so the union-find saw every lighter alternative"
            " first and the cut falls on a component's heaviest links"
        )
    return (
        "edge order: the input's own (--sort off, the profiler's default), so the cut falls wherever"
        " the ceil(n/k) threshold bit and NOT on a component's heaviest links"
    )


def clamp_note(cells):
    """What `--subtract-overhead` clamped, when it clamped enough to matter.

    A correction that drives a column to zero is not a correction, it is the backend telling you the
    column was mostly clock. On the Linux `rdpmc` backend the subtraction is tens of nanoseconds and
    this never fires; on a `clock_gettime` backend at hundreds of nanoseconds a read, `cp_len + 2`
    reads can exceed the whole makespan at small `d`, and every number downstream — the mean, the
    speedup — is then an artefact of the clamp rather than a measurement. So it is said on the
    figure, not only in the table.
    """
    clamped = sum(cell.clamped for cell in cells)
    if not clamped:
        return ""
    shots = sum(cell.shots for cell in cells) or 1
    share = 100.0 * clamped / (2 * shots)
    note = f"{clamped:,} corrected values fell below zero and were clamped to it ({share:.2f}% of columns)"
    if share >= CLAMP_ALARM_PERCENT:
        note = (
            "WARNING: " + note + " — at this rate the correction exceeds the measurement, so the"
            " corrected means and every speedup drawn from them are artefacts of the clamp. Read the"
            " raw figures on this backend."
        )
    return note


def escalation_note(cell):
    return (
        f"{cell.escalated:,} of {cell.shots:,} shots escalated ({100.0 * cell.escalation_rate:.3f}%)"
        " and are NOT in the sparse_k series: an escalating shot stops at the job that truncated, so"
        " it has no makespan. fallback and system cover every shot."
    )


def timer_note(run, args):
    note = f"timer {run.timer}"
    if not run.thread_scoped:
        note += " — WALL CLOCK, timings include off-CPU time"
    if args.subtract_overhead:
        note += (
            f"  ·  overhead-corrected: {run.overhead_ns:.3f} ns x (cp_len + {PRE_INTERIOR_READS})"
            f" off sparse_k per shot and x {INTERVALS_FALLBACK} off fallback; cp_* means are raw"
        )
    else:
        note += "  ·  raw uncorrected ticks (--subtract-overhead applies the per-shot correction)"
    return note


def out_stem(args, run, stem):
    """A figure's filename stem, made unique across runs when they share one output directory.

    Without this, two runs of the same grid written to one `--out-dir` produce the same names and the
    second silently overwrites the first — the failure mode where you compare two campaigns and are
    shown one of them twice. The default (no `--out-dir`, each run writing into its own directory)
    cannot collide, so the prefix is added only when it can.
    """
    if args.out_dir and len(args.paths) > 1:
        return f"{run.name}_{stem}"
    return stem


def save(fig, args, run, stem, theme):
    written = []
    out_dir = args.out_dir or run.path
    os.makedirs(out_dir, exist_ok=True)
    for extension in args.formats:
        out_path = os.path.join(out_dir, f"{out_stem(args, run, stem)}.{extension}")
        fig.savefig(out_path, dpi=args.dpi, facecolor=theme["surface"])
        written.append(out_path)
    plt.close(fig)
    return written


CAPTION_FONT_SIZE = 8.5


# One character of caption text, as a fraction of the font size. matplotlib's default face is DejaVu
# Sans, whose advance over mixed-case prose with spaces averages close to this; `caption_width` is
# what turns it into a wrap column, and the margin below it absorbs the error either way.
CAPTION_CHAR_EM = 0.52


def caption_width(fig, margin_in=0.32):
    """How many characters of caption fit across `fig`.

    `wrap_caption`'s fixed 135 was tuned on the wide `k`-sweeps. A run with a single `k` draws its
    small multiples four inches across, and there the last third of every caption wrapped past the
    right edge of the canvas, where it was cropped without a warning — the figure looked finished and
    was missing the sentence that said what it was. So the column comes from the figure's own width.
    """
    usable_points = max(1.0, (fig.get_figwidth() - margin_in) * 72.0)
    return max(40, int(usable_points / (CAPTION_FONT_SIZE * CAPTION_CHAR_EM)))


def caption_block(fig, theme, caption, top=0.94):
    """One caption under a figure, and the `rect` that leaves exactly room for it.

    The reserve is computed from the figure's own height rather than taken from a constant, because
    these figures range from five inches tall to ten and a fraction that suits one leaves a band of
    dead surface under the other. A line of 8.5pt text with normal leading is about 0.16 inch, so the
    band is `(lines * 0.16 + padding) / height` — capped, so a runaway caption cannot squeeze the
    panels to nothing.
    """
    caption = wrap_caption(caption, width=caption_width(fig))
    lines = caption.count("\n") + 1
    height = fig.get_figheight()
    reserve = min(0.40, (lines * 0.16 + 0.18) / height)
    fig.text(0.012, 0.012, caption, color=theme["text_secondary"], fontsize=CAPTION_FONT_SIZE, va="bottom")
    fig.tight_layout(rect=(0, reserve, 1, top))


PAPER_LABELS = ("pymatching on $\\mathcal{G}$ (fallback)", "SpecMatching, $k$ cores")


def draw_cell_paper(run, cell, edges, args, theme):
    """One cell at single-column width for a paper: the histograms and their means, nothing else.

    No title and no caption block. Everything the long caption states — the shot count, the
    escalations left out of the sparse series, the timer and the k-core model — belongs in the
    paper's own caption, where it can be said once and in the paper's notation.
    """
    fig, ax = plt.subplots(figsize=(3.5, 2.1))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])
    labels = (PAPER_LABELS[0], PAPER_LABELS[1].replace("$k$", f"$k={cell.k:g}$"))
    draw_histograms(
        ax,
        cell,
        edges,
        theme,
        log_x=args.log_x,
        series_labels=labels,
        show_system=False,
        fontsize=7,
        linewidth=1.2,
        unit=r"$\mu$s",
    )
    style_axes(ax, theme, args.log_x, args.log_y)
    ax.tick_params(labelsize=7)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel(r"per-shot latency ($\mu$s)", color=theme["text_secondary"], fontsize=8)
    ax.set_ylabel("shots", color=theme["text_secondary"], fontsize=8)
    # Thousands of shots, so the tick labels stay narrow at column width.
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value / 1000:g}k" if value else "0"))
    legend = ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=2, fontsize=7)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])
    fig.tight_layout(pad=0.2)
    return save(fig, args, run, f"flb_{cell.slug}", theme)


def draw_cell(run, cell, edges, args, theme):
    """One `(d, p, T, k)` cell: the two distributions, their means, and what the caption must state."""
    if args.paper:
        written = draw_cell_paper(run, cell, edges, args, theme)
        fallback, sparse, system = stats_of(cell.fallback), stats_of(cell.sparse), stats_of(cell.system)
        return written, {
            "cell": cell,
            "fallback": fallback,
            "sparse": sparse,
            "system": system,
            "speedup": fallback["mean"] / system["mean"] if system["mean"] else float("nan"),
            "clipped": 0,
        }
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    draw_histograms(ax, cell, edges, theme, log_x=args.log_x)
    style_axes(ax, theme, args.log_x, args.log_y)
    ax.set_xlim(edges[0], edges[-1])
    ax.set_xlabel("per-shot latency (microseconds)", color=theme["text_secondary"], fontsize=10)
    ax.set_ylabel("shots", color=theme["text_secondary"], fontsize=10)
    ax.set_title(cell.title, color=theme["text_primary"], fontsize=13, loc="left", pad=14)
    legend = ax.legend(frameon=False, loc="upper right", fontsize=10)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    fallback = stats_of(cell.fallback)
    sparse = stats_of(cell.sparse)
    system = stats_of(cell.system)
    speedup = fallback["mean"] / system["mean"] if system["mean"] else float("nan")
    clipped = int(
        np.count_nonzero(cell.fallback > edges[-1]) + np.count_nonzero(cell.sparse > edges[-1])
    )
    pieces = cell.histogram_mean("n_pieces_hist")
    depth = cell.histogram_mean("tree_depth_hist")
    caption = (
        f"{fallback['n']:,} shots  ·  speedup {speedup:.2f}x = fallback mean / system mean  ·  "
        f"max {fallback['max']:,.2f} vs {sparse['max']:,.2f} us  ·  "
        f"mean {pieces:.2f} pieces, tree depth {depth:.2f}, cp_len {cell.cp_len:.2f} jobs  ·  "
        f"{clipped:,} beyond the right edge\n"
        f"{escalation_note(cell)}\n"
        "system = min(fallback, sparse_k) on an ordinary shot and fallback alone on an escalating"
        " one, since that decode was already running; it is a reference rule rather than a third"
        f" series because it is a per-shot min of the two that are drawn\n{timer_note(run, args)}"
    )
    caveat = model_caveat(run)
    if caveat:
        caption += f"\n{caveat}"
    clamped = clamp_note([cell])
    if clamped:
        caption += f"\n{clamped}"
    if cell.system_mismatch:
        caption += (
            f"\nWARNING: {cell.system_mismatch:,} rows where the logged ticks_system disagrees with"
            " the definition recomputed from ticks_fallback, ticks_sparse_k and escalated"
        )
    caption_block(fig, theme, caption, top=1.0)
    return save(fig, args, run, f"flb_{cell.slug}", theme), {
        "cell": cell,
        "fallback": fallback,
        "sparse": sparse,
        "system": system,
        "speedup": speedup,
        "clipped": clipped,
    }


def draw_k_sweep(run, key, cells, edges, args, theme):
    """The `k` sweep as small multiples: one panel per `k`, shared bins, shared axes.

    Small multiples rather than `|k|` overlaid curves because the palette will not carry more than
    two hues — see the module docstring. It is also the better read: the question is what the
    `sparse_k` body does as `k` grows, and a panel per `k` answers it without asking anyone to
    untangle four curves.

    `fallback` is repeated in every panel deliberately. It is the same shots decoded the same way at
    every `k`, so it is the control: it should sit still, and a panel where it has moved is telling
    you about the machine rather than about the cut.
    """
    if len(cells) < 2:
        return None
    columns = len(cells)
    fig, axes = plt.subplots(1, columns, figsize=(min(4.6 * columns, 18.0), 5.2), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    fig.patch.set_facecolor(theme["surface"])

    for ax, cell in zip(axes, cells):
        ax.set_facecolor(theme["surface"])
        draw_histograms(ax, cell, edges, theme, log_x=args.log_x, label=False)
        style_axes(ax, theme, args.log_x, args.log_y)
        ax.set_xlim(edges[0], edges[-1])
        speedup = mean_of(cell.fallback) / mean_of(cell.system) if mean_of(cell.system) else float("nan")
        ax.set_title(
            f"k = {cell.k:g}\nsparse_k mean {mean_of(cell.sparse):,.2f} us  ·  speedup {speedup:.2f}x",
            color=theme["text_primary"],
            fontsize=10,
            loc="left",
            pad=10,
        )
        ax.set_xlabel("per-shot latency (microseconds)", color=theme["text_secondary"], fontsize=9)
    axes[0].set_ylabel("shots", color=theme["text_secondary"], fontsize=10)

    fig.suptitle(
        f"{group_title(key)}  —  what k buys",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        y=0.985,
        ha="left",
    )
    # At the figure level rather than inside a panel: the identity is the same in every panel, and a
    # legend box inside one of them sits on that panel's data and on no other's.
    handles, labels = axes[0].get_legend_handles_labels()
    legend = fig.legend(handles, labels, frameon=False, fontsize=9, loc="upper right", bbox_to_anchor=(0.995, 1.0))
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])
    caption = (
        "One panel per k, shared bins and shared axes. The dashed rules are the two means; the dotted"
        " rule is the system mean. fallback is the control — the same shots decoded the same way in"
        " every panel, so it should not move with k.\n"
        + escalation_note(cells[0]).replace(f"{cells[0].escalated:,} of ", "at k = "
                                            f"{cells[0].k:g}, {cells[0].escalated:,} of ")
        + f"\n{timer_note(run, args)}"
    )
    clamped = clamp_note(cells)
    if clamped:
        caption += f"\n{clamped}"
    caveat = model_caveat(run)
    if caveat:
        caption += f"\n{caveat}"
    caption_block(fig, theme, caption)
    return save(fig, args, run, f"flb_ksweep_{group_slug(key)}", theme)


def axis_max(rows):
    """The largest finite value across a grid of bar heights, or 1.0 if there is none.

    Exists because a missing `agg.json` key reads as NaN all the way through to `set_xlim`, where
    matplotlib raises and takes the whole run down — which is how a schema change turns into a
    traceback three figures into a campaign instead of a blank bar on one panel.
    """
    finite = [value for row in rows for value in row if math.isfinite(value)]
    return max(finite) if finite else 1.0


def draw_bar_panel(ax, theme, labels, values, unit="ns"):
    """One horizontal breakdown: one bar per row, one hue, a direct label on each.

    One hue for every bar because the categories are identity, not magnitude, and the row label is
    what carries the identity. A ramp here would double-encode length as lightness. Linear and not
    log, because a bar encodes length from a zero baseline and a log axis would grow a sliver into
    something that looks comparable with the work.
    """
    positions = np.arange(len(labels))
    ax.barh(positions, values, height=0.62, color=theme["series"][1], zorder=2)
    ax.set_yticks(positions)
    ax.set_yticklabels(labels)
    # The row label is the identity, and it is drawn once per row on the leftmost panel; a tick mark
    # on a panel that carries no label is ink with nothing to point at.
    ax.tick_params(axis="y", length=0)
    # Five ticks, not matplotlib's default: at large `d` the labels run to six digits and the default
    # spacing butts them into one another.
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    style_axes(ax, theme, grid_axis="x")
    # A direct label on every bar, which is the one place that is right: a handful of bars, orders of
    # magnitude between them, and the smallest is often the point.
    for position, value in zip(positions, values):
        if not math.isfinite(value):
            continue
        ax.annotate(
            f"{value:,.0f}",
            xy=(value, position),
            ha="left",
            va="center",
            fontsize=8.5,
            color=theme["text_primary"],
            xytext=(4, 0),
            textcoords="offset points",
            zorder=4,
        )
    del unit


def draw_critical_path(run, key, cells, args, theme):
    """Where the makespan goes, and what the manager's share of it is made of.

    Two rows of panels, one column per `k`. The top row is the six critical-path buckets, which sum
    to the `sparse_k` mean by construction. The bottom row is the manager core's own work: the three
    preprocessing regions, which sum to `pre`, plus the dispatches and the combine.

    The `cp_*` means are raw even under `--subtract-overhead` — the split of `cp_len` across the six
    buckets is not logged and inventing one would be worse than saying so. The module means below
    them are corrected at one interval each, so the two rows are captioned separately.
    """
    usable = [cell for cell in cells if cell.agg and cell.agg.get("critical_path_means")]
    if not usable:
        return None
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    columns = len(usable)
    fig, axes = plt.subplots(
        2, columns, figsize=(min(5.0 * columns, 18.0), 8.2), sharex="row", sharey="row", squeeze=False
    )
    fig.patch.set_facecolor(theme["surface"])

    path_values = [[cell.path_mean(name) for name in CRITICAL_PATH] for cell in usable]
    module_values = [[cell.module(name, overhead) for name in MANAGER_MODULES] for cell in usable]
    # Over the finite values only. A bucket this run does not carry comes back as NaN, and a NaN axis
    # limit is a traceback rather than a figure: the reader's job is to draw what is there and leave
    # a visible gap where something is not, not to die because the schema moved.
    path_max = axis_max(path_values)
    module_max = axis_max(module_values)

    for column, cell in enumerate(usable):
        top = axes[0][column]
        top.set_facecolor(theme["surface"])
        draw_bar_panel(top, theme, CRITICAL_PATH, path_values[column])
        top.set_xlim(0, path_max * 1.34)
        total = sum(value for value in path_values[column] if math.isfinite(value))
        sparse_mean = float((cell.agg or {}).get("sparse_k", {}).get("mean", float("nan")))
        agree = math.isfinite(sparse_mean) and abs(total - sparse_mean) <= 1.0
        top.set_title(
            f"k = {cell.k:g}\nsum {total:,.0f} ns"
            + (" = sparse_k mean" if agree else f"  <-- sparse_k mean is {sparse_mean:,.0f}"),
            color=theme["text_primary"],
            fontsize=10,
            loc="left",
            pad=10,
        )

        bottom = axes[1][column]
        bottom.set_facecolor(theme["surface"])
        draw_bar_panel(bottom, theme, MANAGER_MODULES, module_values[column])
        bottom.set_xlim(0, module_max * 1.34)
        pre_sum = sum(cell.module(name, overhead) for name in PRE_MODULES)
        pre = cell.module("pre", overhead)
        pre_agree = math.isfinite(pre) and abs(pre_sum - pre) <= 1.0
        bottom.set_title(
            f"manager core, k = {cell.k:g}\nuf+tree+scatter {pre_sum:,.0f} ns"
            + (" = pre mean" if pre_agree else f"  <-- pre mean is {pre:,.0f}"),
            color=theme["text_primary"],
            fontsize=10,
            loc="left",
            pad=10,
        )
        bottom.set_xlabel("mean ticks (nanoseconds)", color=theme["text_secondary"], fontsize=9)

    # Once per row, not once per panel: the panels share a y axis, so inverting inside the loop
    # inverts the same axis `|k|` times and an even `|k|` silently leaves the rows bottom-up.
    axes[0][0].invert_yaxis()
    axes[1][0].invert_yaxis()

    fig.suptitle(
        f"{group_title(key)}  —  where the makespan goes",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        y=0.99,
        ha="left",
    )
    caption = (
        "Top row: the six buckets of the critical path, which sum to the sparse_k mean by"
        " construction — the path is walked back through whatever determined each start, so every"
        " link is tight and idle is 0. A non-zero idle means the schedule and its own critical path"
        " disagree, which is a bug and not a measurement. build and leaf are the two halves of one"
        " LEAF job: the piece's graph build, then its solve.\n"
        "Bottom row: the manager core's own work. The first three are the preprocessing regions and"
        " sum to pre; dispatch_total is every job assignment and combine is the final reduction."
        " These are serial — they do not shrink with k, and dispatch_total grows with it. There is no"
        " bucket region, and no sort is timed in either edge order.\n"
        f"{edge_order_note(run)}.\n"
        "Means are over NON-ESCALATING shots only, as agg.json accumulates them.\n"
        + (
            f"cp_* are RAW even under --subtract-overhead (the split of cp_len across the six"
            f" buckets is not logged); the module row is corrected at {INTERVALS_MODULE} interval"
            " each.\n"
            if args.subtract_overhead
            else ""
        )
        + timer_note(run, args)
    )
    caption_block(fig, theme, caption, top=0.95)
    return save(fig, args, run, f"flb_path_{group_slug(key)}", theme)


def draw_work(run, key, cells, args, theme):
    """Solver work summed over the `k` cores against the part of it on the critical path.

    This is the question the binary exists to answer. Cutting an oversized component does not make
    the solver do less work — it makes the work divisible — so the payoff is visible only as the gap
    between the two bars, and only if the left bar stays put while the right one falls. A `k` where
    both rise is a `k` where the cut has started costing more than it saves.

    Two hues, which is the budget, and the categories are direct-labelled as well as legended.
    """
    usable = [cell for cell in cells if cell.agg and cell.agg.get("module_means")]
    if len(usable) < 1:
        return None
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    fig, ax = plt.subplots(figsize=(max(6.0, 1.9 * len(usable) + 3.0), 5.2))
    fig.patch.set_facecolor(theme["surface"])
    ax.set_facecolor(theme["surface"])

    positions = np.arange(len(usable), dtype=np.float64)
    total = [sum(cell.module(name, overhead) for name in WORK_TOTAL) for cell in usable]
    on_path = [sum(cell.path_mean(name) for name in WORK_PATH) for cell in usable]
    width = 0.32
    # A surface gap between the two bars of a pair, so the pair reads as a pair of marks rather than
    # as one divided block.
    gap = 0.02
    left = positions - width / 2 - gap
    right = positions + width / 2 + gap
    ax.bar(
        left,
        total,
        width=width,
        color=theme["series"][0],
        label="build + solve + fuse, summed over the k cores",
        zorder=2,
    )
    ax.bar(
        right,
        on_path,
        width=width,
        color=theme["series"][1],
        label="build + solve + fuse, on the critical path",
        zorder=2,
    )
    for position, values in ((left, total), (right, on_path)):
        for at, value in zip(position, values):
            if not math.isfinite(value):
                continue
            ax.annotate(
                f"{value:,.0f}",
                xy=(at, value),
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=theme["text_primary"],
                xytext=(0, 3),
                textcoords="offset points",
                zorder=4,
            )
    ax.set_xticks(positions)
    ax.set_xticklabels([f"k = {cell.k:g}" for cell in usable])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
    style_axes(ax, theme, log_y=args.log_y)
    # Headroom for two things at once: the direct label above the tallest bar, and the legend box in
    # the corner above it. Without it the tallest bar grows straight through both.
    tallest = max((value for value in total + on_path if math.isfinite(value)), default=1.0)
    if not args.log_y:
        ax.set_ylim(0, tallest * 1.34)
    ax.set_ylabel("mean ticks (nanoseconds)", color=theme["text_secondary"], fontsize=10)
    ax.set_title(
        f"{group_title(key)}  —  work done against work on the path",
        color=theme["text_primary"],
        fontsize=13,
        loc="left",
        pad=14,
    )
    legend = ax.legend(frameon=False, loc="upper right", fontsize=9)
    for text in legend.get_texts():
        text.set_color(theme["text_primary"])

    ratios = "  ·  ".join(
        f"k={cell.k:g}: {(t / p):.2f}x" if p else f"k={cell.k:g}: n/a"
        for cell, t, p in zip(usable, total, on_path)
    )
    caption = (
        "Left bar: build_sum + solve_sum + fuse_sum, the solver time the shot costs in total. Right"
        " bar: cp_build + cp_leaf + cp_fuse, the part of it that lands on the critical path. Cutting"
        " does not reduce the work, it makes the work divisible, so the payoff is the gap — and only"
        " if the left bar stays put while the right one falls.\n"
        "The graph build is in both bars. It is per-piece work on a solver core, so it divides with"
        " everything else; a left bar that omitted it would be measuring a cut that costs nothing to"
        " set up.\n"
        f"work / path: {ratios}\n"
        "Neither bar includes the manager's serial share; see the critical-path figure for that."
        " Means over NON-ESCALATING shots only. cp_* are raw even under --subtract-overhead.\n"
        f"{timer_note(run, args)}"
    )
    caveat = model_caveat(run)
    if caveat:
        caption += f"\n{caveat}"
    caption_block(fig, theme, caption, top=1.0)
    return save(fig, args, run, f"flb_work_{group_slug(key)}", theme)


def draw_structure(run, key, cells, args, theme):
    """The cut itself: how many pieces, how deep the tree, how wide a fusion, how many releases.

    Rows are the four statistics, columns are `k`, which is the small-multiple form the two-hue
    budget forces and also the right one — the question is how each distribution *moves* with `k`,
    and a grid answers it by letting the eye run along a row.

    Counts come from `agg.json`, where the profiler histogrammed them over the whole cell.
    """
    usable = [cell for cell in cells if cell.agg]
    if not usable:
        return None
    columns = len(usable)
    fig, axes = plt.subplots(
        len(STRUCTURE), columns, figsize=(min(4.2 * columns, 18.0), 10.0), sharey="row", squeeze=False
    )
    fig.patch.set_facecolor(theme["surface"])
    drew = False

    for row, (name, label) in enumerate(STRUCTURE):
        # Shared x per row so the eye can run along it; the four statistics have unrelated supports,
        # so sharing x down a column would be meaningless.
        highs = []
        for cell in usable:
            values, counts = cell.histogram(name)
            if values.size:
                highs.append(float(values.max()))
        high = max(highs) if highs else 1.0
        for column, cell in enumerate(usable):
            ax = axes[row][column]
            ax.set_facecolor(theme["surface"])
            values, counts = cell.histogram(name)
            if values.size:
                ax.bar(values, counts, width=0.85, color=theme["series"][1], zorder=2)
                drew = True
            ax.set_xlim(-0.6, high + 0.6)
            ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=6))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
            style_axes(ax, theme, log_y=args.log_y)
            if row == 0:
                ax.set_title(f"k = {cell.k:g}", color=theme["text_primary"], fontsize=10, loc="left", pad=8)
            if column == 0:
                ax.set_ylabel(label, color=theme["text_secondary"], fontsize=9)
            if row == len(STRUCTURE) - 1:
                ax.set_xlabel("value", color=theme["text_secondary"], fontsize=9)
            mean = cell.histogram_mean(name)
            if math.isfinite(mean):
                ax.annotate(
                    f"mean {mean:.2f}",
                    xy=(0.97, 0.92),
                    xycoords="axes fraction",
                    ha="right",
                    va="top",
                    fontsize=8.5,
                    color=theme["text_primary"],
                    zorder=4,
                )
    if not drew:
        plt.close(fig)
        return None

    fig.suptitle(
        f"{group_title(key)}  —  the cut, and the tree it produces",
        color=theme["text_primary"],
        fontsize=13,
        x=0.012,
        y=0.99,
        ha="left",
    )
    caption = (
        "Counts over the cell's measured shots, from agg.json. Pieces are what the bounded union-find"
        " left after refusing any edge that would take a piece past ceil(n/k); the tree is the"
        " cheapest-boundary-first greedy over them; a fusion's boundary is how many crossing edges it"
        " turns interior; a release is a region whose dummy boundary was raised out from under it.\n"
        f"{edge_order_note(run)} — every panel here is downstream of that choice.\n"
        "One hue: these are counts of one thing per panel, so identity is carried by the row label"
        " and not by colour. A shot with no defects contributes a 0 to every row."
    )
    caption_block(fig, theme, caption, top=0.95)
    return save(fig, args, run, f"flb_structure_{group_slug(key)}", theme)


def fmt(value, width, digits=2):
    return f"{value:>{width},.{digits}f}" if isinstance(value, float) and math.isfinite(value) else f"{'-':>{width}}"


def table_lines(run, results, args):
    """Everything the figures show, as text, so the plots are never the only read.

    One list of lines, so the same table can go to the terminal and to the file beside the figures
    without the two drifting apart.
    """
    lines = []
    lines.append(f"fusion_lb_profiler — {run.path}")
    lines.append(f"  git {run.meta.get('git_hash', '?')}  ·  stim {run.meta.get('stim_version', '?')}")
    lines.append(
        f"  timer {run.timer}  ·  thread-scoped {'yes' if run.thread_scoped else 'NO — WALL CLOCK'}"
        f"  ·  overhead {run.overhead_ns:.3f} ns/read"
    )
    lines.append(f"  {edge_order_note(run)}")

    def first_word(key):
        return run.meta.get(key, "?").split()[0] if run.meta.get(key) else "?"

    lines.append(
        f"  seed {first_word('seed')}  ·  alpha {run.alpha:g}"
        f"  ·  shots/cell {first_word('shots_per_cell')}"
        f"  ·  warmup/cell {first_word('warmup_shots_per_cell')}"
        f"  ·  cells {len(results)}"
    )
    lines.append(
        "  correction: "
        + (
            f"{run.overhead_ns:.3f} ns x (cp_len + {PRE_INTERIOR_READS}) off sparse_k per shot,"
            f" x {INTERVALS_FALLBACK} off fallback, x {INTERVALS_MODULE} off each module mean;"
            " cp_* left raw"
            if args.subtract_overhead
            else "none — raw uncorrected ticks"
        )
    )
    for key in ("sparse_k", "system", "horizon_definition"):
        note = run.meta.get(key)
        if note:
            lines.append(f"  {key} = {note}")
    caveat = model_caveat(run)
    if caveat:
        for part in caveat.split("  ·  "):
            lines.append(f"  ! {part}")
    lines.append(
        "  ! an escalating shot has no makespan (ticks_sparse_k = 0) and is excluded from every"
        " sparse_k and cp_* figure; fallback and system cover every shot"
    )
    lines.append("")

    # ---- latency.
    header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'shots':>9} {'escal':>7} {'escal%':>7} "
        f"{'fallback':>10} {'sparse_k':>10} {'system':>10} {'speedup':>8} "
        f"{'fb max':>10} {'sp max':>10} {'sys p99':>10} {'sys p99.9':>10}"
    )
    lines.append(header)
    lines.append("-" * len(header))
    for result in results:
        cell = result["cell"]
        system_agg = (cell.agg or {}).get("system", {})
        p99 = system_agg.get("p99")
        p999 = system_agg.get("p99_9")
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} {cell.shots:>9,} "
            f"{cell.escalated:>7,} {100.0 * cell.escalation_rate:>7.3f} "
            f"{fmt(result['fallback']['mean'], 10)} {fmt(result['sparse']['mean'], 10)} "
            f"{fmt(result['system']['mean'], 10)} {fmt(result['speedup'], 8)} "
            f"{fmt(result['fallback']['max'], 10)} {fmt(result['sparse']['max'], 10)} "
            f"{fmt(p99 / 1000.0 if p99 is not None else float('nan'), 10)} "
            f"{fmt(p999 / 1000.0 if p999 is not None else float('nan'), 10)}"
        )
    lines.append("")
    lines.append("  all latencies in microseconds; percentile columns are agg.json's, over RAW ticks")
    lines.append("  sparse_k columns are over non-escalating shots only; fallback and system over all")

    # ---- the critical path, and the identity it has to satisfy.
    lines.append("")
    cp_header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} "
        + " ".join(f"{'cp_' + name:>12}" for name in CRITICAL_PATH)
        + f" {'sum':>12} {'sparse_k':>12} {'cp_len':>7}"
    )
    lines.append(cp_header)
    lines.append("-" * len(cp_header))
    for result in results:
        cell = result["cell"]
        if not (cell.agg or {}).get("critical_path_means"):
            continue
        values = [cell.path_mean(name) for name in CRITICAL_PATH]
        total = sum(value for value in values if math.isfinite(value))
        sparse_mean = float((cell.agg or {}).get("sparse_k", {}).get("mean", float("nan")))
        flag = "" if not math.isfinite(sparse_mean) or abs(total - sparse_mean) <= 1.0 else "  <-- MISMATCH"
        idle = cell.path_mean("idle")
        if math.isfinite(idle) and idle != 0.0:
            flag += "  <-- cp_idle NOT ZERO"
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
            + " ".join(fmt(value, 12, 1) for value in values)
            + f" {fmt(total, 12, 1)} {fmt(sparse_mean, 12, 1)} {fmt(cell.cp_len, 7)}{flag}"
        )
    lines.append("")
    lines.append("  critical-path means in RAW nanoseconds, from agg.json, over non-escalating shots")
    lines.append("  the six sum to the sparse_k mean by construction; MISMATCH means they no longer do")
    lines.append("  cp_build and cp_leaf are the two halves of one LEAF job: the build, then the solve")
    lines.append("  cp_idle is 0 by construction — every link of the walked-back path is tight")

    # ---- the manager core.
    lines.append("")
    mm_header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} "
        + " ".join(f"{name:>13}" for name in MANAGER_MODULES)
        + f" {'pre sum':>13} {'pre':>13}"
    )
    lines.append(mm_header)
    lines.append("-" * len(mm_header))
    overhead = run.overhead_ns if args.subtract_overhead else 0.0
    for result in results:
        cell = result["cell"]
        if not (cell.agg or {}).get("module_means"):
            continue
        values = [cell.module(name, overhead) for name in MANAGER_MODULES]
        pre_sum = sum(cell.module(name, overhead) for name in PRE_MODULES)
        pre = cell.module("pre", overhead)
        flag = "" if not math.isfinite(pre) or abs(pre_sum - pre) <= 1.0 else "  <-- MISMATCH"
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
            + " ".join(fmt(value, 13, 1) for value in values)
            + f" {fmt(pre_sum, 13, 1)} {fmt(pre, 13, 1)}{flag}"
        )
    lines.append("")
    lines.append("  manager-core means in nanoseconds; uf+tree+scatter = pre by construction")
    lines.append("  no bucket region, and no sort is timed in either edge order; see the edge order above")

    # ---- work against path. One triple of columns per kind of solver work, so that the build —
    # which used to be the manager's whole-instance rebuild inside scatter — is visible as its own
    # line rather than folded into the solve.
    lines.append("")
    work_header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'build_sum':>12} {'build_max':>12} {'cp_build':>12} "
        f"{'solve_sum':>12} {'solve_max':>12} {'cp_leaf':>12} "
        f"{'fuse_sum':>12} {'fuse_max':>12} {'cp_fuse':>12} {'extract_max':>12} {'work/path':>10}"
    )
    lines.append(work_header)
    lines.append("-" * len(work_header))
    for result in results:
        cell = result["cell"]
        if not (cell.agg or {}).get("module_means"):
            continue
        total = sum(cell.module(name, overhead) for name in WORK_TOTAL)
        on_path = sum(cell.path_mean(name) for name in WORK_PATH)
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
            f"{fmt(cell.module('build_sum', overhead), 12, 1)} {fmt(cell.module('build_max', overhead), 12, 1)} "
            f"{fmt(cell.path_mean('build'), 12, 1)} "
            f"{fmt(cell.module('solve_sum', overhead), 12, 1)} {fmt(cell.module('solve_max', overhead), 12, 1)} "
            f"{fmt(cell.path_mean('leaf'), 12, 1)} "
            f"{fmt(cell.module('fuse_sum', overhead), 12, 1)} {fmt(cell.module('fuse_max', overhead), 12, 1)} "
            f"{fmt(cell.path_mean('fuse'), 12, 1)} {fmt(cell.module('extract_max', overhead), 12, 1)} "
            f"{fmt(total / on_path if on_path else float('nan'), 10)}"
        )
    lines.append("")
    lines.append("  build/solve/fuse_* are the k cores' work; cp_* is the part of it on the path")
    lines.append("  work/path = (build_sum+solve_sum+fuse_sum) / (cp_build+cp_leaf+cp_fuse)")

    # ---- the cut's structure.
    lines.append("")
    st_header = (
        f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'pieces':>9} {'depth':>9} {'boundary':>9} {'released':>9}"
    )
    lines.append(st_header)
    lines.append("-" * len(st_header))
    for result in results:
        cell = result["cell"]
        lines.append(
            f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
            + " ".join(fmt(cell.histogram_mean(name), 9) for name, _ in STRUCTURE)
        )
    lines.append("")
    lines.append("  means of agg.json's structure histograms: pieces, fusion-tree depth,")
    lines.append("  largest fusion boundary (crossing edges turned interior), regions released")

    # ---- what the run asserted about itself.
    gates = [
        result
        for result in results
        if (result["cell"].agg or {}).get("verify_shots") is not None
    ]
    if gates:
        lines.append("")
        gate_header = (
            f"{'d':>4} {'p':>8} {'T':>5} {'k':>3} {'verified':>9} {'div_fused':>10} {'(leaf':>6}"
            f"{'/fuse)':>7} {'div_mono':>9} {'determ':>7} {'thresh':>7} {'cut':>5} {'feas':>5}"
            f" {'mgrwr':>6} {'bindep':>7} {'bfused':>7}"
        )
        lines.append(gate_header)
        lines.append("-" * len(gate_header))
        for result in gates:
            agg = result["cell"].agg
            cell = result["cell"]
            determinism = {1: "pass", 0: "FAIL"}.get(agg.get("determinism_check"), "-")
            lines.append(
                f"{cell.d:>4.0f} {cell.p:>8g} {cell.T:>5g} {cell.k:>3.0f} "
                f"{agg.get('verify_shots', 0):>9,} "
                f"{agg.get('escalation_divergence_fused_only', 0):>10,} "
                f"{agg.get('escalation_divergence_fused_only_in_leaf', 0):>6,}"
                f"{agg.get('escalation_divergence_fused_only_in_fusion', 0):>7,} "
                f"{agg.get('escalation_divergence_monolithic_only', 0):>9,} "
                f"{determinism:>7} "
                f"{agg.get('threshold_violations', 0):>7,} "
                f"{agg.get('cut_on_heaviest_violations', 0):>5,} "
                f"{agg.get('feasibility_violations', 0):>5,} "
                f"{agg.get('manager_graph_write_violations', 0):>6,} "
                f"{agg.get('build_independence_violations', 0):>7,} "
                f"{agg.get('built_before_fused_violations', 0):>7,}"
            )
        lines.append("")
        lines.append("  verified = shots cross-checked against the monolithic solve of the whole of H")
        lines.append("  div_* = shots where only one of the two escalated; a RATE difference, not an")
        lines.append("          exactness one — both directions are safe, see the profiler's deviation 9")
        lines.append("  thresh/cut/feas = the profiler's own §11 invariants; anything but 0 is a bug")
        if run.sorted_edges is False:
            lines.append("  cut = 0 here because it was NOT CHECKED: the two cut-on-heaviest tests are claims")
            lines.append("          about a sorted input, and this run did not sort. thresh and feas did run")
        lines.append("  mgrwr/bindep/bfused = the per-piece build's own invariants: no solver-node write")
        lines.append("          on the manager between tree and scatter, leaves in reverse order leaving")
        lines.append("          byte-identical state, and both children built before a fusion. They are")
        lines.append("          debug-build checks, so a release run reports 0 because none ran")

    mismatches = sum(result["cell"].system_mismatch for result in results)
    if mismatches:
        lines.append("")
        lines.append(
            f"  WARNING: {mismatches:,} rows where the logged ticks_system disagrees with the"
            " definition recomputed from ticks_fallback, ticks_sparse_k and escalated"
        )
    clamped = clamp_note([result["cell"] for result in results])
    if clamped:
        lines.append("")
        lines.append(f"  {clamped}")
    return lines


def warn_about_wall_clock(run):
    if run.thread_scoped:
        return
    print(
        f"\n  WARNING: {run.path} was measured with a WALL CLOCK, backend"
        f" '{run.timer}'.\n  Every tick includes time the thread spent off the CPU and is an upper"
        " bound only.\n",
        file=sys.stderr,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument(
        "paths",
        nargs="+",
        help="fusion_lb_profiler --out directories (shots_d*_p*_T*_k*.csv + agg.json + run.log)",
    )
    parser.add_argument("--out-dir", help="where the figures go; default is inside each run directory")
    parser.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--bins", type=int, default=90)
    parser.add_argument("--log-x", action="store_true", help="log latency axis; the escalation tail is long")
    parser.add_argument("--log-y", action="store_true", help="log count axis, which is where a thin tail shows up")
    parser.add_argument(
        "--x-max-percentile",
        type=float,
        default=99.5,
        help="where the latency axis stops; a view control, and the caption counts the shots past it",
    )
    parser.add_argument(
        "--subtract-overhead",
        action="store_true",
        help=f"subtract timer_overhead_ns x (cp_len + {PRE_INTERIOR_READS}) from sparse_k per shot"
        " and x 1 from fallback; off by default, and the caption always says which is on the page",
    )
    parser.add_argument("--d", help="only these distances, comma separated")
    parser.add_argument("--p", help="only these error rates, comma separated")
    parser.add_argument("--T", help="only these horizons, comma separated")
    parser.add_argument("--k", help="only these core counts, comma separated")
    parser.add_argument("--no-cells", action="store_true", help="skip the per-cell histograms")
    parser.add_argument("--no-ksweep", action="store_true", help="skip the k-sweep small multiples")
    parser.add_argument("--no-path", action="store_true", help="skip the critical-path breakdown")
    parser.add_argument("--no-work", action="store_true", help="skip the work-against-path figure")
    parser.add_argument("--no-structure", action="store_true", help="skip the cut-structure grid")
    parser.add_argument(
        "--paper",
        action="store_true",
        help="per-cell histograms at single-column width with no title or caption block, for a paper"
        " whose own caption carries the context; combine with --formats pdf",
    )
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--table", action="store_true", help="print the summary table to stdout as well")
    args = parser.parse_args()

    args.formats = [item.strip() for item in args.formats.split(",") if item.strip()]
    args.d = parse_number_list(args.d, float)
    args.p = parse_number_list(args.p, float)
    args.T = parse_number_list(args.T, float)
    args.k = parse_number_list(args.k, float)
    theme = THEMES[args.theme]

    failed = False
    for path in args.paths:
        try:
            run = load(path)
        except (ValueError, OSError, json.JSONDecodeError) as error:
            print(f"  {path}: {error}", file=sys.stderr)
            failed = True
            continue
        warn_about_wall_clock(run)
        print(f"{run.path}  ({len(run.cell_files)} cell files)")
        results = []
        for key, cells in stream_groups(run, args):
            if not cells:
                continue
            x_max, edges = axis_limits(cells, args)
            for cell in cells:
                if not args.no_cells:
                    written, result = draw_cell(run, cell, edges, args, theme)
                    for out_path in written:
                        print(f"  wrote {out_path}")
                else:
                    result = {
                        "cell": cell,
                        "fallback": stats_of(cell.fallback),
                        "sparse": stats_of(cell.sparse),
                        "system": stats_of(cell.system),
                        "clipped": 0,
                    }
                    result["speedup"] = (
                        result["fallback"]["mean"] / result["system"]["mean"]
                        if result["system"]["mean"]
                        else float("nan")
                    )
                results.append(result)
            for draw, skip in (
                (lambda: draw_k_sweep(run, key, cells, edges, args, theme), args.no_ksweep),
                (lambda: draw_critical_path(run, key, cells, args, theme), args.no_path),
                (lambda: draw_work(run, key, cells, args, theme), args.no_work),
                (lambda: draw_structure(run, key, cells, args, theme), args.no_structure),
            ):
                if skip:
                    continue
                for out_path in draw() or ():
                    print(f"  wrote {out_path}")
            # The group's arrays are released here, before the next one is read.
            for cell in cells:
                cell.release()

        if not results:
            print("  no cells matched the selectors; nothing to plot", file=sys.stderr)
            failed = True
            continue
        lines = table_lines(run, results, args)
        out_dir = args.out_dir or run.path
        os.makedirs(out_dir, exist_ok=True)
        table_path = os.path.join(out_dir, f"{out_stem(args, run, 'fusion_lb_latency')}.txt")
        with open(table_path, "w") as handle:
            handle.write("\n".join(lines) + "\n")
        print(f"  wrote {table_path}")
        if args.table:
            print()
            print("\n".join(lines))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
