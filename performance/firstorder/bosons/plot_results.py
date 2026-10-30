#!/usr/bin/env python3
"""firstorder/bosons/plot_results.py —— 一阶冻结核 vs 顶点移位改善标度图。"""
import csv, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
SPEC = [
    ("mixedtime.csv", "mixedtime.png", "dt",
     f"Mixed time (N_tau = N_t): frozen-kernel schemes vs strict Trotter\n{'Fermionic' if 'fermions' in HERE else 'Bosonic'} Toulouse,  G>(t),  first-order frozen kernel  vs  strict Trotter 1/2"),
    ("realtime.csv", "realtime.png", "dt",
     "Real time (bosonic Toulouse, A(t)):  frozen-kernel schemes vs strict Trotter\nfirst-order assembly  |  vertex-shifted improvement  |  strict order 1/2"),
    ("imagtime.csv", "imagtime.png", "dtau",
     "Imaginary time (bosonic Toulouse, G(tau)):  frozen-kernel schemes vs strict Trotter\nfirst-order ring  |  frozen-cell improvement  |  strict order 1/2"),
]
COLS = {"err_frozen1": ("o-", "#d62728", 1.8, 7, "frozen kernel, 1st-order assembly"),
        "err_shifted": ("s-", "#ff7f0e", 1.8, 7, "vertex-shifted improvement (real-time axis)"),
        "err_strict1": ("^:", "#7f7f7f", 1.4, 6, "strict chain, order 1 (reference)"),
        "err_strict2": ("s--", "#1f77b4", 1.5, 6, "strict chain, order 2 (reference)")}

for csv_name, png_name, xlabel, title in SPEC:
    path = os.path.join(RESULTS, csv_name)
    if not os.path.exists(path):
        print(f"skip {csv_name}"); continue
    rows = list(csv.DictReader(open(path)))
    step = [float(r["dt"]) if "dt" in r else float(r["dtau"]) for r in rows]
    fig, ax = plt.subplots(figsize=(6.6, 4.9))
    for col, (style, color, lw, ms, label) in COLS.items():
        if col not in rows[0]:
            continue
        vals = [float(r[col]) for r in rows]
        if max(vals) < 1e-13:
            continue
        ax.loglog(step, vals, style, color=color, lw=lw, ms=ms, label=label)
    s0, e0 = step[0], [float(r["err_frozen1"]) for r in rows][0]
    span = (max(step)/min(step))**0.45
    for p, a in ((1, 0.45), (2, 0.45)):
        xs = [min(step)/span, max(step)*span]
        ax.loglog(xs, [e0*(x/s0)**p for x in xs], "k-", lw=0.7, alpha=a)
        ax.annotate(rf"$\propto h^{p}$", xy=(xs[1], e0*(xs[1]/s0)**p),
                    xytext=(6, 2), textcoords="offset points", fontsize=11, alpha=0.75)
    ax.set_xlabel(rf"${xlabel}$", fontsize=12); ax.set_ylabel("max error vs exact", fontsize=12)
    ax.set_title(title, fontsize=10)
    ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6); ax.legend(fontsize=8.5)
    ax.invert_xaxis(); fig.tight_layout()
    out = os.path.join(RESULTS, png_name)
    fig.savefig(out, dpi=170); plt.close(fig); print(f"wrote {out}")
