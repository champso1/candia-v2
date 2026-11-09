#!/usr/bin/env python3
"""Plots data/CT18NLO.dat and data/CT18NNLO.dat from subtraction-pdfs as subtraction-pdfs.pdf,
one page per figure of the first three appendix pages of the S-ACOT-MPS working note:
bottom PDFs at fixed scale, versus scale, and near threshold."""

import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpecFromSubplotSpec

SETS = [("CT18NLO", "CT18 NLO"), ("CT18NNLO", "CT18 NNLO")]
# column number of each quantity in the data file
QUANTITIES = {
    "fb":     (2, r"$f_b$",                       "#1f5fa6", "-"),
    "ft_lo":  (3, r"$\widetilde{f}_b^{\rm LO}$",  "#c8442b", (0, (5, 2.5))),
    "ft_nlo": (4, r"$\widetilde{f}_b^{\rm NLO}$", "#16856f", (0, (1, 1.5))),
    "df_lo":  (5, r"$\delta f_b^{\rm LO}$",       "#b8620f", (0, (6, 2, 1.5, 2))),
    "df_nlo": (6, r"$\delta f_b^{\rm NLO}$",      "#6b3fc4", (0, (6, 1.5, 1.5, 1.5, 1.5, 1.5))),
}
MAIN = ["fb", "ft_lo", "ft_nlo"]
RESIDUALS = ["df_lo", "df_nlo"]
FLOOR = 1e-8  # positive-only log axis of the residual panels

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 9,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.3,
    "legend.fontsize": 8,
    "legend.frameon": True,
    "legend.fancybox": False,
    "legend.edgecolor": "0.3",
})


def load(name):
    data = np.loadtxt(f"data/{name}.dat")
    with open(f"data/{name}.dat") as f:
        mb = float(f.readlines()[1].split("m_b =")[1].split()[0])
    return data, mb


def at_scale(data, q):
    rows = data[np.isclose(data[:, 0], q, rtol=1e-9)]
    return rows[np.argsort(rows[:, 1])]


def at_x(data, x):
    rows = data[np.isclose(data[:, 1], x, rtol=1e-9)]
    return rows[np.argsort(rows[:, 0])]


def curve(ax, abscissa, rows, key, floor=None, label=True):
    col, tex, color, style = QUANTITIES[key]
    y = rows[:, col]
    ax.plot(abscissa, y, color=color, linestyle=style, label=tex if label else None)


def style(ax):
    ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=8)


def header(fig, title, subtitle):
    height = fig.get_figheight()  # place the lines a fixed distance (in inches) from the top
    fig.text(0.08, 1.0 - 0.25/height, title, fontsize=13, ha="left", va="top")
    fig.text(0.08, 1.0 - 0.52/height, subtitle, fontsize=8.5, ha="left", va="top")


def fixed_scale(pdf, sets):
    fig = plt.figure(figsize=(8.5, 8.0))
    outer = fig.add_gridspec(2, 2, left=0.10, right=0.97, top=0.88, bottom=0.10, hspace=0.32, wspace=0.18)
    for row, q in enumerate([10.0, 100.0]):
        for col, (name, label) in enumerate(SETS):
            data, _ = sets[name]
            rows = at_scale(data, q)
            x = rows[:, 1]
            inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[row, col], height_ratios=[2, 1], hspace=0.0)
            top = fig.add_subplot(inner[0])
            bottom = fig.add_subplot(inner[1], sharex=top)
            for key in MAIN:
                curve(top, x, rows, key)
            for key in RESIDUALS:
                curve(bottom, x, rows, key, floor=FLOOR)
            top.set(xscale="log", yscale="log", xlim=(1e-5, 1), ylim=(1e-5, 10))
            bottom.set(yscale="log", ylim=(FLOOR, 1), xlabel=r"$x$")
            top.set_title(rf"{label} $\cdot$ $Q = {q:g}\,$GeV", loc="left", fontsize=10)
            plt.setp(top.get_xticklabels(), visible=False)
            if col == 0:
                top.set_ylabel(r"$xf_b,\ x\widetilde{f}_b$")
                bottom.set_ylabel(r"$x\,\delta f_b$")
            if row == 0 and col == 0:
                top.legend(loc="lower left")
                bottom.legend(loc="upper right")
            for ax in (top, bottom):
                style(ax)
    header(fig, "Bottom PDFs at fixed scale",
           r"Candia-v2 evolution of CT18 NLO and NNLO $\cdot$ Paper Fig. 4 kinematics")
    pdf.savefig(fig)
    plt.close(fig)


def versus_scale(pdf, sets):
    fig, axes = plt.subplots(2, 3, figsize=(8.5, 7.0), gridspec_kw=dict(left=0.09, right=0.97, top=0.87, bottom=0.10, hspace=0.38, wspace=0.28))
    for row, (name, label) in enumerate(SETS):
        data, mb = sets[name]
        for col, xval in enumerate([1e-4, 1e-2, 0.1]):
            ax = axes[row, col]
            rows = at_x(data, xval)
            q = rows[:, 0]
            ax.axhline(0.0, color="0.6", linewidth=0.6)
            ax.axvline(mb, color="0.4", linewidth=0.6, linestyle=":")
            for key in MAIN + RESIDUALS:
                curve(ax, q, rows, key)
            ax.set(xscale="log", xlim=(4.5, 1e4), xlabel=r"$Q$ [GeV]")
            xtex = r"10^{-4}" if xval == 1e-4 else (r"10^{-2}" if xval == 1e-2 else "0.1")
            ax.set_title(rf"{label} $\cdot$ $x = {xtex}$", loc="left", fontsize=10)
            if col == 0:
                ax.set_ylabel(r"$xf_b,\ x\widetilde{f}_b,\ x\,\delta f_b$")
            if row == 0 and col == 0:
                ax.legend(loc="upper left")
            style(ax)
    header(fig, "Bottom PDFs versus scale",
           r"Candia-v2 evolution of CT18 NLO and NNLO $\cdot$ Paper Fig. 5 kinematics")
    pdf.savefig(fig)
    plt.close(fig)


def near_threshold(pdf, sets):
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.6), sharey=True, gridspec_kw=dict(left=0.09, right=0.97, top=0.80, bottom=0.17, wspace=0.12))
    for col, (name, label) in enumerate(SETS):
        data, _ = sets[name]
        rows = at_scale(data, 5.0)
        ax = axes[col]
        ax.axhline(0.0, color="0.6", linewidth=0.6)
        for key in MAIN + RESIDUALS:
            curve(ax, rows[:, 1], rows, key)
        ax.set(xscale="log", xlim=(1e-5, 1), xlabel=r"$x$")
        ax.set_title(label, loc="left", fontsize=10)
        if col == 0:
            ax.set_ylabel(r"$xf_b,\ x\widetilde{f}_b,\ x\,\delta f_b$")
            ax.legend(loc="upper right")
        style(ax)
    header(fig, "Bottom PDFs near threshold",
           r"Candia-v2 evolution of CT18 NLO and NNLO $\cdot$ $Q = 5\,$GeV $\cdot$ Paper Fig. 6 kinematics")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    sets = {name: load(name) for name, _ in SETS}
    with PdfPages("subtraction-pdfs.pdf") as pdf:
        fixed_scale(pdf, sets)
        versus_scale(pdf, sets)
        near_threshold(pdf, sets)
    print("wrote subtraction-pdfs.pdf")


if __name__ == "__main__":
    main()
