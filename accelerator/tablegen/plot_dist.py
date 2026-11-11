#!/usr/bin/env python3
"""CCDFs of the per-shot latency distributions written by sweep_dist.py.

Usage:  conda run -n syndrilla python plot_dist.py <dist-dir>

Reads <dist-dir>/dist_d<d>_p<p>.csv and writes <dist-dir>/latency_dist.png (300 dpi) and
latency_dist.pdf: one panel per (d, p), one curve per lane count, P(cycles > x) on log-log axes.

Needs matplotlib. Words are set in Gill Sans, numbers in cmss10 through mathtext; cmss10 ships
with matplotlib and its x-height is within 1.3% of Gill Sans at 8 pt, so no size matching is
needed. If Gill Sans is missing, matplotlib falls back to its default sans.
"""

import csv
import os
import sys

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullLocator, StrMethodFormatter

DISTANCES = [17, 25, 31]
NOISES = [0.0005, 0.001, 0.003]
LANES = [16, 32, 64, 128]

# Okabe-Ito, colourblind-safe; the dash pattern repeats the distinction so colour is never alone.
STYLE = {
    16: ("#0072B2", "-"),
    32: ("#E69F00", (0, (5, 1.5))),
    64: ("#009E73", (0, (1.5, 1.2))),
    128: ("#CC79A7", (0, (4, 1.2, 1, 1.2))),
}
TEXT_PT = 8


def set_style():
    mpl.rcParams.update(
        {
            "font.family": "Gill Sans",
            "font.size": TEXT_PT,
            "mathtext.fontset": "custom",
            "mathtext.rm": "cmss10",
            "mathtext.it": "cmss10",
            "mathtext.bf": "cmss10",
            "mathtext.cal": "cmss10",
            "xtick.labelsize": TEXT_PT,
            "ytick.labelsize": TEXT_PT,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "legend.frameon": False,
        }
    )


def ccdf(values):
    """x and P(X > x), as a right-continuous step: the k-th largest sits at k / n."""
    vals = sorted(values)
    n = len(vals)
    return vals, [1.0 - (i + 1) / n for i in range(n)]


def read_point(dist_dir, d, p):
    path = os.path.join(dist_dir, "dist_d%d_p%g.csv" % (d, p))
    with open(path) as f:
        rows = [{k: int(v) for k, v in r.items()} for r in csv.DictReader(f)]
    if not rows:
        sys.exit("%s: no rows" % path)
    return rows


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    dist_dir = sys.argv[1]
    set_style()

    fig, axes = plt.subplots(3, 3, figsize=(7.0, 6.4))
    n_shots = None
    for row, d in enumerate(DISTANCES):
        for col, p in enumerate(NOISES):
            ax = axes[row][col]
            rows = read_point(dist_dir, d, p)
            n_shots = len(rows)
            for lanes in LANES:
                colour, dash = STYLE[lanes]
                x, y = ccdf(r["cycles_L%d" % lanes] for r in rows)
                ax.step(x, y, where="post", color=colour, linestyle=dash, linewidth=1.1)
            ax.set_xscale("log")
            ax.set_yscale("log")
            ax.set_ylim(1.0 / n_shots, 1.4)
            ax.set_title("d = $%d$,  p = $%g$" % (d, p), pad=4)
            ax.grid(True, which="major", linewidth=0.3, alpha=0.4)
            ax.tick_params(length=2.5)
            # A panel can span less than a decade, where the default log locator falls back to
            # linear-spaced ticks that collide. Pinning x to 1/2/5 per decade keeps every panel
            # readable and in the digit face; minor ticks would only reintroduce the collision.
            ax.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0), numticks=99))
            ax.xaxis.set_minor_locator(NullLocator())
            for axis in (ax.xaxis, ax.yaxis):
                axis.set_major_formatter(StrMethodFormatter("${x:g}$"))
            if col == 0:
                ax.set_ylabel("P(cycles > x)")
            if row == len(DISTANCES) - 1:
                ax.set_xlabel("H-build latency (cycles)")

    handles = [
        Line2D(
            [],
            [],
            color=STYLE[L][0],
            linestyle=STYLE[L][1],
            linewidth=1.1,
            label="$%d$ lanes" % L,
        )
        for L in LANES
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=len(LANES),
        bbox_to_anchor=(0.5, 1.0),
        columnspacing=2.0,
        handlelength=2.6,
    )
    fig.suptitle(
        "h_construct latency distribution, legacy ball layout, per-shot capacity",
        y=1.045,
    )
    fig.tight_layout()

    png = os.path.join(dist_dir, "latency_dist.png")
    pdf = os.path.join(dist_dir, "latency_dist.pdf")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    print("wrote %s and %s (n = %d shots per panel)" % (png, pdf, n_shots))


if __name__ == "__main__":
    main()
