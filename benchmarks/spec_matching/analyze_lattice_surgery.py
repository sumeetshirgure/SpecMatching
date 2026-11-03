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

"""Tables for what `lattice_surgery_profiler` logs (docs/lattice_surgery_profiler.md §2.4).

    python benchmarks/spec_matching/analyze_lattice_surgery.py RUN_DIR [RUN_DIR ...] [--out DIR]

A run directory holds `agg.json`, `run.log` and, per cell, `shots_*.csv`, `jobs_*.csv` and
`phase2_*.csv`. Several directories may be given — the sweep script writes one per corpus — and
their cells are pooled into one table. A cell is identified by the `# N=..,k=..,...` comment line
of its CSVs, never by parsing the filename.

What is reported, per cell, following §2.4:

 1. **Seam-fusion depth** — the fraction of shots with a spanning chain, and mean `seam_depth`.
 2. **Latency** — `t_sys` against B1 (`t_strip_pm`, stock on the whole strip) and B2
    (`t_global_escalate`, the paper's §4.5 semantics), mean and p99.9, per shot and per round
    (divided by `d`, the number of rounds). B3 (`t_cube_pm_max`) is a latency reference only.
 3. **Fallback statistics** — the strip-level truncation rate and the per-cube rate, against the
    `q` of an `N = 1` cell at the same `(d, p, T, k, orientation)` when the run has one (that cell is
    this generator's `d × d × d` memory experiment), and the independent-cube prediction
    `1 - (1 - q)^N`.
 4. **Locality** — over truncated shots: final `Ω` in detectors and in cubes, `(t1)`/`(t2)`
    triggers, reset components and defects, hand-offs, and how often `Ω` touched more than three
    cubes.
 5. **Exactness** — `--verify`'s weight and observable mismatch counts, and phase-2 conflicts.

Rows are raw, uncorrected thread-CPU nanoseconds; `run.log` has `timer_overhead_ns`. Nothing here
subtracts it: a makespan is a sum over however many jobs land on the critical path, and inventing a
count would be worse than saying so.
"""

import argparse
import json
import math
import os
import sys
from glob import glob

import numpy as np
import pandas as pd

KEY = ("N", "k", "d", "p", "T", "W", "orientation")
T_COLUMNS = (
    "t_sys_ns",
    "t_spec_ns",
    "t_strip_pm_ns",
    "t_global_escalate_ns",
    "t_cube_pm_max_ns",
    "t_fb_max_ns",
    "t_p2_ns",
)
CP_COLUMNS = (
    "cp_accel",
    "cp_mgr",
    "cp_leaf",
    "cp_fuse",
    "cp_seam",
    "cp_extract",
    "cp_gm",
    "cp_fb",
    "cp_p2",
    "cp_link",
    "cp_idle",
)
SHOT_COLUMNS = T_COLUMNS + CP_COLUMNS + (
    "n_def",
    "n_cubes_truncated",
    "n_spanning",
    "seam_depth",
    "omega_nodes",
    "omega_cubes",
    "n_t1",
    "n_t2",
    "n_closure_steps",
    "n_reset_components",
    "n_reset_defects",
    "n_handoff",
    "handoff_nodes",
    "n_moved_defects",
    "p2_conflict",
    "weight",
    "ref_weight",
    "obs_xor",
    "ref_obs",
)


def parse_header(path):
    """The `# N=..,k=..` comment line of a cell's CSV, as a dict."""
    with open(path) as f:
        line = f.readline().strip()
    if not line.startswith("#"):
        return None
    out = {}
    for item in line[1:].split(","):
        if "=" in item:
            key, value = item.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def cell_key(meta):
    return (
        int(meta["N"]),
        int(meta["k"]),
        int(meta["d"]),
        float(meta["p"]),
        float(meta["T"]),
        float(meta["W"]),
        meta["orientation"],
    )


def load_cells(run_dirs):
    """Every cell of every run directory: its metadata, its shot rows and its `agg.json` entry."""
    cells = {}
    for run_dir in run_dirs:
        agg_entries = {}
        agg_path = os.path.join(run_dir, "agg.json")
        if os.path.exists(agg_path):
            with open(agg_path) as f:
                for entry in json.load(f).get("cells", []):
                    key = (entry["N"], entry["k"], entry["d"], float(entry["p"]), float(entry["T"]),
                           float(entry["W"]), entry["orientation"])
                    agg_entries[key] = entry
        for path in sorted(glob(os.path.join(run_dir, "shots_*.csv"))):
            meta = parse_header(path)
            if meta is None:
                print(f"skipping {path}: no header comment", file=sys.stderr)
                continue
            key = cell_key(meta)
            shots = pd.read_csv(path, comment="#", usecols=lambda c: c in SHOT_COLUMNS)
            if key in cells:
                print(f"warning: cell {key} appears twice; keeping {path}", file=sys.stderr)
            cells[key] = {
                "meta": meta,
                "shots": shots,
                "agg": agg_entries.get(key, {}),
                "path": path,
                "phase2": path.replace("shots_", "phase2_"),
                "jobs": path.replace("shots_", "jobs_"),
                "run_dir": run_dir,
            }
    return cells


def p999(values):
    values = np.asarray(values)
    if values.size == 0:
        return float("nan")
    return float(np.quantile(values, 0.999, method="lower"))


def mean(values):
    values = np.asarray(values, dtype=float)
    return float(values.mean()) if values.size else float("nan")


def summarise_cell(key, cell, q_reference):
    s = cell["shots"]
    N, k, d, p, T, W, orientation = key
    truncated = s["n_cubes_truncated"] > 0
    nt = s[~truncated]
    tr = s[truncated]
    row = dict(zip(KEY, key))
    row["shots"] = len(s)
    row["shots_truncated"] = int(truncated.sum())
    row["strip_trunc_rate"] = truncated.mean()
    row["cube_trunc_rate"] = s["n_cubes_truncated"].sum() / (len(s) * N)
    q = q_reference.get((k, d, p, T, orientation))
    row["q_ref"] = q if q is not None else float("nan")
    row["strip_rate_independent"] = 1 - (1 - q) ** N if q is not None else float("nan")
    row["spanning_frac"] = (s["n_spanning"] > 0).mean()
    row["seam_depth_mean"] = s["seam_depth"].mean()
    row["seam_depth_ge2_frac"] = (s["seam_depth"] >= 2).mean()
    for col in ("t_sys_ns", "t_strip_pm_ns", "t_global_escalate_ns", "t_cube_pm_max_ns"):
        name = col[:-3]
        row[f"{name}_mean_us"] = mean(s[col]) / 1e3
        row[f"{name}_p999_us"] = p999(s[col]) / 1e3
        row[f"{name}_mean_us_per_round"] = row[f"{name}_mean_us"] / d
    row["t_spec_mean_us"] = mean(nt["t_spec_ns"]) / 1e3
    row["t_spec_p999_us"] = p999(nt["t_spec_ns"]) / 1e3
    row["speedup_vs_B1"] = row["t_strip_pm_mean_us"] / row["t_sys_mean_us"]
    row["speedup_vs_B2"] = row["t_global_escalate_mean_us"] / row["t_sys_mean_us"]
    row["t_fb_max_mean_us"] = mean(s["t_fb_max_ns"]) / 1e3
    # Truncated shots only: what local resolution cost, and how local it stayed.
    row["trunc_t_sys_mean_us"] = mean(tr["t_sys_ns"]) / 1e3
    row["trunc_t_strip_pm_mean_us"] = mean(tr["t_strip_pm_ns"]) / 1e3
    row["trunc_t_p2_mean_us"] = mean(tr["t_p2_ns"]) / 1e3
    row["trunc_t_fb_max_mean_us"] = mean(tr["t_fb_max_ns"]) / 1e3
    row["omega_nodes_mean"] = mean(tr["omega_nodes"])
    row["omega_cubes_mean"] = mean(tr["omega_cubes"])
    row["omega_cubes_max"] = int(tr["omega_cubes"].max()) if len(tr) else 0
    row["omega_cubes_gt3"] = int((s["omega_cubes"] > 3).sum())
    row["triggers_mean"] = mean(tr["n_t1"] + tr["n_t2"])
    row["t1_mean"] = mean(tr["n_t1"])
    row["t2_mean"] = mean(tr["n_t2"])
    row["closure_steps_mean"] = mean(tr["n_closure_steps"])
    row["reset_components_mean"] = mean(tr["n_reset_components"])
    row["reset_defects_mean"] = mean(tr["n_reset_defects"])
    row["handoff_shots"] = int((s["n_handoff"] > 0).sum())
    row["handoff_nodes_mean"] = mean(tr["handoff_nodes"])
    row["moved_defects_mean"] = mean(s["n_moved_defects"])
    row["p2_conflicts"] = int(s["p2_conflict"].sum())
    row["weight_mismatches"] = int((s["weight"] != s["ref_weight"]).sum())
    row["obs_mismatches"] = int((s["obs_xor"] != s["ref_obs"]).sum())
    row["verified"] = bool(cell["agg"].get("verify_shots", 0))
    for col in CP_COLUMNS:
        row[f"{col}_mean_ns"] = mean(nt[col])
    return row


def q_references(cells):
    """`q` per `(k, d, p, T, orientation)` from the `N = 1` cells, the strip's memory reference."""
    out = {}
    for key, cell in cells.items():
        N, k, d, p, T, W, orientation = key
        if N != 1:
            continue
        s = cell["shots"]
        out[(k, d, p, T, orientation)] = float((s["n_cubes_truncated"] > 0).mean())
    # A k with no N=1 cell borrows any k's: truncation does not depend on k, only the schedule does.
    for key in list(out):
        k, d, p, T, orientation = key
        for other in cells:
            if other[2:5] == (d, p, T) and other[6] == orientation and (other[1], d, p, T, orientation) not in out:
                out[(other[1], d, p, T, orientation)] = out[key]
    return out


def summary_frame(cells):
    q = q_references(cells)
    rows = [summarise_cell(key, cells[key], q) for key in sorted(cells)]
    return pd.DataFrame(rows)


def fmt(value, digits=2):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "–"
    if isinstance(value, (int, np.integer)):
        return f"{value:,}"
    if abs(value) != 0 and (abs(value) < 1e-2 or abs(value) >= 1e5):
        return f"{value:.2e}"
    return f"{value:,.{digits}f}"


def markdown_table(df, columns, headers):
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("--:" for _ in headers) + "|"]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join(fmt(r[c]) if not isinstance(r[c], str) else r[c] for c in columns) + " |")
    return "\n".join(lines)


def write_markdown(df, path, run_dirs):
    out = ["# lattice_surgery_profiler summary", ""]
    out.append("Run directories: " + ", ".join(f"`{r}`" for r in run_dirs))
    out.append("")
    out.append("Raw thread-CPU ticks, uncorrected. `t_spec` is over non-truncated shots; every other latency is over "
               "all shots. Per round = per shot / d.")
    ident = ["N", "k", "d", "p", "T", "W", "orientation"]
    sections = [
        ("Latency (§2.4.2)", ident + ["shots", "t_sys_mean_us", "t_sys_p999_us", "t_sys_mean_us_per_round",
                                     "t_strip_pm_mean_us", "t_global_escalate_mean_us", "t_cube_pm_max_mean_us",
                                     "speedup_vs_B1", "speedup_vs_B2"],
         ident + ["shots", "t_sys µs", "t_sys p99.9 µs", "t_sys µs/round", "B1 µs", "B2 µs", "B3 µs",
                  "vs B1", "vs B2"]),
        ("Truncation and seams (§2.4.1, §2.4.3)", ident + ["strip_trunc_rate", "cube_trunc_rate", "q_ref",
                                                           "strip_rate_independent", "spanning_frac",
                                                           "seam_depth_ge2_frac"],
         ident + ["strip rate", "per-cube rate", "q (N=1)", "1-(1-q)^N", "spanning shots", "depth ≥ 2"]),
        ("Truncated shots: cost and locality (§2.4.4)", ident + ["shots_truncated", "trunc_t_sys_mean_us",
                                                                "trunc_t_strip_pm_mean_us", "trunc_t_fb_max_mean_us",
                                                                "trunc_t_p2_mean_us", "omega_nodes_mean",
                                                                "omega_cubes_mean", "omega_cubes_max", "triggers_mean",
                                                                "reset_defects_mean", "handoff_shots",
                                                                "omega_cubes_gt3"],
         ident + ["truncated", "t_sys µs", "B1 µs", "phase 1 µs", "phase 2 µs", "Ω nodes", "Ω cubes", "Ω cubes max",
                  "t1+t2", "reset defects", "hand-off shots", "Ω > 3 cubes"]),
        ("Exactness (§2.4.5)", ident + ["shots", "weight_mismatches", "obs_mismatches", "p2_conflicts"],
         ident + ["shots", "weight mismatches", "obs mismatches", "phase-2 conflicts"]),
        ("Critical path, non-truncated shots, mean ns", ident + [f"{c}_mean_ns" for c in CP_COLUMNS],
         ident + [c[3:] for c in CP_COLUMNS]),
    ]
    for title, cols, heads in sections:
        out += ["", f"## {title}", "", markdown_table(df, cols, heads)]
    if not df["verified"].all():
        out += ["", "Note: some cells were run without `--verify`; their mismatch columns compare against B1 anyway, "
                    "since B1 runs on every shot."]
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_dirs", nargs="+")
    parser.add_argument("--out", help="where summary.csv and summary.md go; default: the first run directory")
    args = parser.parse_args()
    run_dirs = []
    for r in args.run_dirs:
        # A sweep root holds one sub-directory per corpus; take them all.
        subs = sorted(os.path.dirname(p) for p in glob(os.path.join(r, "*", "agg.json")))
        run_dirs += subs if subs and not os.path.exists(os.path.join(r, "agg.json")) else [r]
    cells = load_cells(run_dirs)
    if not cells:
        sys.exit("no cells found")
    df = summary_frame(cells)
    out = args.out or run_dirs[0]
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "summary.csv"), index=False)
    write_markdown(df, os.path.join(out, "summary.md"), run_dirs)
    print(f"{len(df)} cells -> {out}/summary.csv, {out}/summary.md")
    bad = df[df["weight_mismatches"] > 0]
    if len(bad):
        print("WARNING: weight mismatches against B1 in", len(bad), "cells", file=sys.stderr)


if __name__ == "__main__":
    main()
