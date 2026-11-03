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

"""Figures for what `lattice_surgery_profiler` logs.

    python benchmarks/spec_matching/plot_lattice_surgery.py RUN_DIR [RUN_DIR ...] [--out DIR]
                                                            [--formats png,pdf] [--theme light|dark]

Reads cells through `analyze_lattice_surgery.load_cells`, so a sweep root (one sub-directory per
corpus) or any list of run directories works. Figures, one set per `(orientation, k, T, W)`:

  * `latency_vs_N_*`   — mean `t_sys` per round against B1 (stock on the strip), with B2 (global
                          escalation) as a grey reference rule; one panel per `(d, p)`. The design's
                          headline comparison (§8.6). A `_p999` twin draws the 99.9th percentiles.
  * `speedup_vs_N_*`   — `B1 / t_sys` and `B2 / t_sys`. The B2 gap is what local resolution buys and
                          should grow with `N`.
  * `truncation_vs_N_*` — the strip-level truncation rate against `1 - (1 - q)^N`, the
                          independent-cube prediction from the `N = 1` cell (§2.4.3). Log scale.
  * `spanning_vs_N_*`  — the fraction of shots whose seam fusion needs depth ≥ 2 (§2.4.1).
  * `critical_path_*`  — where the non-truncated makespan goes, one panel per `N`, one bar per
                          component, single hue.

and per cell:

  * `hist_*`     — `t_sys` against B1 per shot.
  * `locality_*` — over truncated shots: cubes `Ω` touched, `(t1) + (t2)` triggers, reset defects,
                    and `|Ω|` (§2.4.4). Only drawn for cells with a truncated shot.

Colour: the corpus palette of `plot_fusion_lb.py`, whose two series slots were checked with the
dataviz validator against both themes' surfaces. No panel carries more than those two hues; a third
quantity is a grey reference rule, and anything with more categories is a small multiple. Identity
is never colour alone: every series is in a legend and the line ends are direct-labelled.
"""

import argparse
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_lattice_surgery import CP_COLUMNS, load_cells, summary_frame  # noqa: E402

THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text_primary": "#0b0b0b",
        "text_secondary": "#52514e",
        "grid": "#dcdcd8",
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

SYS_LABEL = "t_sys: speculative + local resolution"
B1_LABEL = "B1: stock pymatching on the strip"
B2_LABEL = "B2: global escalation"


def style_axes(ax, theme, log_y=False, grid_axis="y"):
    if log_y:
        ax.set_yscale("log")
    ax.set_facecolor(theme["surface"])
    ax.grid(axis=grid_axis, color=theme["grid"], linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(theme["grid"])
    ax.tick_params(colors=theme["text_secondary"], labelsize=9)


def new_figure(theme, rows, cols, width=3.6, height=2.8):
    fig, axes = plt.subplots(rows, cols, figsize=(width * cols, height * rows), squeeze=False)
    fig.patch.set_facecolor(theme["surface"])
    return fig, axes


def finish(fig, theme, title, path_base, formats, dpi):
    fig.suptitle(title, color=theme["text_primary"], fontsize=11, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    for fmt in formats:
        fig.savefig(f"{path_base}.{fmt}", dpi=dpi, facecolor=theme["surface"])
    plt.close(fig)


def panel_title(ax, theme, text):
    ax.set_title(text, color=theme["text_primary"], fontsize=10, loc="left")


def axis_labels(ax, theme, x, y):
    ax.set_xlabel(x, color=theme["text_secondary"], fontsize=9)
    ax.set_ylabel(y, color=theme["text_secondary"], fontsize=9)


def end_label(ax, theme, x, y, text):
    """A direct label at a line's last point, in secondary ink (identity is the line beside it)."""
    ax.annotate(text, (x, y), xytext=(4, 0), textcoords="offset points", va="center", fontsize=8,
                color=theme["text_secondary"])


def group_keys(df):
    return sorted({(r.orientation, r.k, r.T, r.W) for r in df.itertuples()})


def tag(orientation, k, T, W):
    return f"{orientation}_k{k}_T{T:g}_W{W:g}"


def dp_panels(sub):
    return sorted({(r.d, r.p) for r in sub.itertuples()})


def plot_latency(df, theme, out, formats, dpi, stat):
    col = "mean_us" if stat == "mean" else "p999_us"
    for orientation, k, T, W in group_keys(df):
        sub = df[(df.orientation == orientation) & (df.k == k) & (df["T"] == T) & (df.W == W)]
        panels = dp_panels(sub)
        if not panels:
            continue
        fig, axes = new_figure(theme, 1, len(panels))
        for ax, (d, p) in zip(axes[0], panels):
            cell = sub[(sub.d == d) & (sub.p == p)].sort_values("N")
            style_axes(ax, theme)
            x = cell.N.to_numpy()
            ys = cell[f"t_sys_{col}"].to_numpy() / d
            yb1 = cell[f"t_strip_pm_{col}"].to_numpy() / d
            yb2 = cell[f"t_global_escalate_{col}"].to_numpy() / d
            ax.plot(x, yb2, color=theme["text_secondary"], linewidth=1.2, linestyle="--", label=B2_LABEL, zorder=2)
            ax.plot(x, yb1, color=theme["series"][1], linewidth=2, marker="o", markersize=5, label=B1_LABEL, zorder=3)
            ax.plot(x, ys, color=theme["series"][0], linewidth=2, marker="o", markersize=5, label=SYS_LABEL, zorder=4)
            if len(x):
                end_label(ax, theme, x[-1], yb1[-1], f"B1 {yb1[-1]:.2f}")
                end_label(ax, theme, x[-1], ys[-1], f"t_sys {ys[-1]:.2f}")
            ax.set_xscale("log", base=2)
            ax.set_xticks(x)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.set_ylim(bottom=0)
            panel_title(ax, theme, f"d={d}, p={p:g}")
            axis_labels(ax, theme, "patches N", f"{'mean' if stat == 'mean' else 'p99.9'} µs per round")
        axes[0][0].legend(fontsize=8, frameon=False, labelcolor=theme["text_secondary"], loc="upper left")
        name = "latency_vs_N" + ("" if stat == "mean" else "_p999")
        finish(fig, theme, f"Strip latency per round, {orientation}, k={k}, T={T:g}, W={W:g}",
               os.path.join(out, f"{name}_{tag(orientation, k, T, W)}"), formats, dpi)


def plot_speedup(df, theme, out, formats, dpi):
    for orientation, k, T, W in group_keys(df):
        sub = df[(df.orientation == orientation) & (df.k == k) & (df["T"] == T) & (df.W == W)]
        panels = dp_panels(sub)
        fig, axes = new_figure(theme, 1, len(panels))
        for ax, (d, p) in zip(axes[0], panels):
            cell = sub[(sub.d == d) & (sub.p == p)].sort_values("N")
            style_axes(ax, theme)
            x = cell.N.to_numpy()
            ax.axhline(1, color=theme["text_secondary"], linewidth=1, linestyle="--", zorder=1)
            ax.plot(x, cell.speedup_vs_B1, color=theme["series"][1], linewidth=2, marker="o", markersize=5,
                    label="B1 / t_sys", zorder=3)
            ax.plot(x, cell.speedup_vs_B2, color=theme["series"][0], linewidth=2, marker="o", markersize=5,
                    label="B2 / t_sys", zorder=4)
            if len(x):
                end_label(ax, theme, x[-1], cell.speedup_vs_B1.iloc[-1], f"{cell.speedup_vs_B1.iloc[-1]:.2f}×")
                end_label(ax, theme, x[-1], cell.speedup_vs_B2.iloc[-1], f"{cell.speedup_vs_B2.iloc[-1]:.2f}×")
            ax.set_xscale("log", base=2)
            ax.set_xticks(x)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.set_ylim(bottom=0)
            panel_title(ax, theme, f"d={d}, p={p:g}")
            axis_labels(ax, theme, "patches N", "mean speedup (×)")
        axes[0][0].legend(fontsize=8, frameon=False, labelcolor=theme["text_secondary"], loc="upper left")
        finish(fig, theme, f"Speedup of t_sys, {orientation}, k={k}, T={T:g}, W={W:g}",
               os.path.join(out, f"speedup_vs_N_{tag(orientation, k, T, W)}"), formats, dpi)


def plot_truncation(df, theme, out, formats, dpi):
    for orientation, k, T, W in group_keys(df):
        sub = df[(df.orientation == orientation) & (df.k == k) & (df["T"] == T) & (df.W == W)]
        panels = dp_panels(sub)
        fig, axes = new_figure(theme, 1, len(panels))
        for ax, (d, p) in zip(axes[0], panels):
            cell = sub[(sub.d == d) & (sub.p == p)].sort_values("N")
            x = cell.N.to_numpy()
            measured = cell.strip_trunc_rate.replace(0, np.nan)
            predicted = cell.strip_rate_independent.replace(0, np.nan)
            if measured.isna().all() and predicted.isna().all():
                # A rate of zero has no place on a log axis; say so rather than draw an empty frame.
                style_axes(ax, theme)
                ax.text(0.5, 0.5, f"no truncated shots\nin {int(cell.shots.min()):,} per cell", ha="center",
                        va="center", transform=ax.transAxes, color=theme["text_secondary"], fontsize=9)
                ax.set_xticks([])
                ax.set_yticks([])
                panel_title(ax, theme, f"d={d}, p={p:g}")
                continue
            style_axes(ax, theme, log_y=True)
            ax.plot(x, predicted, color=theme["series"][1], linewidth=2, linestyle="--",
                    marker="s", markersize=5, label="1 − (1 − q)^N, q from N = 1", zorder=3)
            ax.plot(x, measured, color=theme["series"][0], linewidth=2, marker="o", markersize=5,
                    label="measured: any cube truncated", zorder=4)
            ax.set_xscale("log", base=2)
            ax.set_xticks(x)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            panel_title(ax, theme, f"d={d}, p={p:g}")
            axis_labels(ax, theme, "patches N", "fraction of shots")
        for ax in axes[0]:
            if ax.get_legend_handles_labels()[0]:
                ax.legend(fontsize=8, frameon=False, labelcolor=theme["text_secondary"], loc="upper left")
                break
        finish(fig, theme, f"Strip-level truncation rate, {orientation}, k={k}, T={T:g}, W={W:g}",
               os.path.join(out, f"truncation_vs_N_{tag(orientation, k, T, W)}"), formats, dpi)


def plot_spanning(df, theme, out, formats, dpi):
    for orientation, k, T, W in group_keys(df):
        sub = df[(df.orientation == orientation) & (df.k == k) & (df["T"] == T) & (df.W == W) & (df.N > 1)]
        panels = dp_panels(sub)
        if not panels:
            continue
        fig, axes = new_figure(theme, 1, len(panels))
        for ax, (d, p) in zip(axes[0], panels):
            cell = sub[(sub.d == d) & (sub.p == p)].sort_values("N")
            style_axes(ax, theme)
            x = cell.N.to_numpy()
            ax.plot(x, cell.seam_depth_ge2_frac, color=theme["series"][0], linewidth=2, marker="o", markersize=5)
            for xi, yi in zip(x, cell.seam_depth_ge2_frac):
                end_label(ax, theme, xi, yi, f"{yi:.3f}")
            ax.set_xscale("log", base=2)
            ax.set_xticks(x)
            ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            ax.set_ylim(0, max(0.05, float(cell.seam_depth_ge2_frac.max()) * 1.2))
            panel_title(ax, theme, f"d={d}, p={p:g}")
            axis_labels(ax, theme, "patches N", "shots with seam depth ≥ 2")
        finish(fig, theme, f"Spanning chains (sequential seam fusion), {orientation}, k={k}, T={T:g}, W={W:g}",
               os.path.join(out, f"spanning_vs_N_{tag(orientation, k, T, W)}"), formats, dpi)


def plot_critical_path(df, theme, out, formats, dpi):
    names = [c[3:] for c in CP_COLUMNS]
    for (orientation, k, T, W) in group_keys(df):
        sub = df[(df.orientation == orientation) & (df.k == k) & (df["T"] == T) & (df.W == W)]
        for d, p in dp_panels(sub):
            cell = sub[(sub.d == d) & (sub.p == p)].sort_values("N")
            fig, axes = new_figure(theme, 1, len(cell), width=3.0, height=3.2)
            vmax = max(float(cell[[f"{c}_mean_ns" for c in CP_COLUMNS]].max().max()) / 1e3, 1e-9)
            for ax, (_, r) in zip(axes[0], cell.iterrows()):
                style_axes(ax, theme, grid_axis="x")
                values = np.array([r[f"{c}_mean_ns"] for c in CP_COLUMNS]) / 1e3
                y = np.arange(len(names))
                ax.barh(y, values, color=theme["series"][0], height=0.6, zorder=2)
                ax.set_yticks(y)
                ax.set_yticklabels(names, fontsize=8)
                ax.invert_yaxis()
                ax.set_xlim(0, vmax * 1.25)
                for yi, v in zip(y, values):
                    if v > 0:
                        ax.annotate(f"{v:.2f}", (v, yi), xytext=(3, 0), textcoords="offset points", va="center",
                                    fontsize=7, color=theme["text_secondary"])
                panel_title(ax, theme, f"N={int(r.N)}  (t_spec {r.t_spec_mean_us:.1f} µs)")
                axis_labels(ax, theme, "mean µs on the critical path", "")
            finish(fig, theme, f"Critical path of non-truncated shots, d={d}, p={p:g}, {orientation}, k={k}, "
                               f"T={T:g}, W={W:g}",
                   os.path.join(out, f"critical_path_d{d}_p{p:g}_{tag(orientation, k, T, W)}"), formats, dpi)


def cell_tag(key):
    N, k, d, p, T, W, orientation = key
    return f"N{N}_k{k}_d{d}_p{p:g}_T{T:g}_W{W:g}_{orientation}"


def plot_histograms(cells, theme, out, formats, dpi, bins):
    for key, cell in sorted(cells.items()):
        s = cell["shots"]
        d = key[2]
        fig, axes = new_figure(theme, 1, 1, width=6.0, height=3.2)
        ax = axes[0][0]
        style_axes(ax, theme, log_y=True)
        a = s["t_sys_ns"].to_numpy() / 1e3
        b = s["t_strip_pm_ns"].to_numpy() / 1e3
        hi = np.quantile(np.concatenate([a, b]), 0.9995)
        edges = np.linspace(0, hi, bins + 1)
        for values, colour, label in ((b, theme["series"][1], B1_LABEL), (a, theme["series"][0], SYS_LABEL)):
            counts, _ = np.histogram(np.clip(values, 0, hi), bins=edges)
            ax.stairs(counts, edges, color=colour, alpha=0.4, fill=True, zorder=2)
            ax.stairs(counts, edges, color=colour, linewidth=2, zorder=3, label=f"{label} (mean {values.mean():.1f})")
            ax.axvline(values.mean(), color=colour, linewidth=1, linestyle="--", zorder=4)
        ax.legend(fontsize=8, frameon=False, labelcolor=theme["text_secondary"])
        axis_labels(ax, theme, f"µs per shot (last bin holds the tail; {d} rounds)", "shots")
        finish(fig, theme, "Per-shot latency, " + cell_tag(key).replace("_", " "),
               os.path.join(out, f"hist_{cell_tag(key)}"), formats, dpi)


def plot_locality(cells, theme, out, formats, dpi):
    for key, cell in sorted(cells.items()):
        s = cell["shots"]
        tr = s[s["n_cubes_truncated"] > 0]
        if tr.empty:
            continue
        quantities = (
            ("omega_cubes", "cubes Ω touches", True),
            (None, "(t1) + (t2) triggers", True),
            ("n_reset_defects", "reset defects", False),
            ("omega_nodes", "Ω, detectors", False),
        )
        fig, axes = new_figure(theme, 1, 4, width=3.0, height=2.8)
        for ax, (col, label, discrete) in zip(axes[0], quantities):
            style_axes(ax, theme)
            values = (tr["n_t1"] + tr["n_t2"]) if col is None else tr[col]
            values = values.to_numpy()
            if discrete:
                vals, counts = np.unique(values, return_counts=True)
                ax.bar(vals, counts, color=theme["series"][0], width=0.6, zorder=2)
                ax.set_xticks(vals)
            else:
                ax.hist(values, bins=min(30, max(5, len(np.unique(values)))), color=theme["series"][0], zorder=2)
            axis_labels(ax, theme, label, "truncated shots")
        finish(fig, theme, f"Locality of phase 2 over {len(tr)} truncated shots, " + cell_tag(key).replace("_", " "),
               os.path.join(out, f"locality_{cell_tag(key)}"), formats, dpi)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("run_dirs", nargs="+")
    parser.add_argument("--out", help="where the figures go; default: <first run dir>/figures")
    parser.add_argument("--formats", default="png", help="comma-separated: png,pdf,svg")
    parser.add_argument("--dpi", type=int, default=160)
    parser.add_argument("--theme", choices=sorted(THEMES), default="light")
    parser.add_argument("--bins", type=int, default=80)
    parser.add_argument("--no-cells", action="store_true", help="skip the per-cell histograms and locality figures")
    args = parser.parse_args()
    from glob import glob

    run_dirs = []
    for r in args.run_dirs:
        subs = sorted(os.path.dirname(p) for p in glob(os.path.join(r, "*", "agg.json")))
        run_dirs += subs if subs and not os.path.exists(os.path.join(r, "agg.json")) else [r]
    cells = load_cells(run_dirs)
    if not cells:
        sys.exit("no cells found")
    df = summary_frame(cells)
    out = args.out or os.path.join(run_dirs[0], "figures")
    os.makedirs(out, exist_ok=True)
    theme = THEMES[args.theme]
    formats = [f.strip() for f in args.formats.split(",") if f.strip()]
    plot_latency(df, theme, out, formats, args.dpi, "mean")
    plot_latency(df, theme, out, formats, args.dpi, "p999")
    plot_speedup(df, theme, out, formats, args.dpi)
    plot_truncation(df, theme, out, formats, args.dpi)
    plot_spanning(df, theme, out, formats, args.dpi)
    plot_critical_path(df, theme, out, formats, args.dpi)
    if not args.no_cells:
        plot_histograms(cells, theme, out, formats, args.dpi, args.bins)
        plot_locality(cells, theme, out, formats, args.dpi)
    print(f"figures -> {out}")


if __name__ == "__main__":
    main()
