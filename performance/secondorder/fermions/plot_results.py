#!/usr/bin/env python3
"""secondorder/fermions/plot_results.py —— 从 results/*.csv 画误差-步长标度图（PNG）。

运行：python3 secondorder/fermions/plot_results.py
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

SPEC = [
    ("realtime.csv", "realtime.png", "dt",
     "Real time: Weyl-K 2nd order (frozen kernel + K side)\nFermionic Toulouse,  A(t),  T = 2.0"),
    ("imagtime.csv", "imagtime.png", "dtau",
     "Imaginary time: Weyl-K 2nd order (frozen kernel + K side)\nFermionic Toulouse,  G(tau),  beta = 4.0"),
    ("mixedtime.csv", "mixedtime.png", "dt",
     "Mixed time (N_tau = N_t): Weyl-K 2nd order (frozen kernel + K side)\nFermionic Toulouse,  G>(t),  beta = 4.0,  T = 2.0"),
]


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    step = [float(r["dt"]) if "dt" in r else float(r["dtau"]) for r in rows]
    return {
        "step": step,
        "weylk": [float(r["err_weylk"]) for r in rows],
        "strict2": [float(r["err_strict2"]) for r in rows],
    }


def main():
    for csv_name, png_name, xlabel, title in SPEC:
        path = os.path.join(RESULTS, csv_name)
        if not os.path.exists(path):
            print(f"skip {csv_name}: not found (run the .jl benchmark first)")
            continue
        d = read_csv(path)
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.loglog(d["step"], d["weylk"], "o-", color="#d62728", lw=1.8, ms=7,
                  label="Weyl-K 2nd order (frozen kernel, K side)")
        ax.loglog(d["step"], d["strict2"], "s--", color="#1f77b4", lw=1.5, ms=6,
                  label="strict chain, order 2 (reference)")
        # 斜率参考线（取首个数据点定标）
        s0, e0 = d["step"][0], d["weylk"][0]
        span = (max(d["step"]) / min(d["step"])) ** 0.45
        ref2 = [e0 * (s / s0) ** 2 for s in (min(d["step"]) / span, max(d["step"]) * span)]
        ref1 = [e0 * (s / s0) ** 1 for s in (min(d["step"]) / span, max(d["step"]) * span)]
        ax.loglog([min(d["step"]) / span, max(d["step"]) * span], ref2, "k-", lw=0.7,
                  alpha=0.45)
        ax.loglog([min(d["step"]) / span, max(d["step"]) * span], ref1, "k-", lw=0.7,
                  alpha=0.45)
        # 在参考线末端标注斜率
        x_end = max(d["step"]) * span
        ax.annotate(r"$\propto h^{2}$", xy=(x_end, e0 * (x_end / s0) ** 2),
                    xytext=(6, 2), textcoords="offset points", fontsize=11, alpha=0.75)
        ax.annotate(r"$\propto h^{1}$", xy=(x_end, e0 * (x_end / s0) ** 1),
                    xytext=(6, -10), textcoords="offset points", fontsize=11, alpha=0.75)
        ax.set_xlabel(rf"${xlabel}$", fontsize=12)
        ax.set_ylabel("max error vs exact", fontsize=12)
        ax.set_title(title, fontsize=10.5)
        ax.grid(True, which="both", ls=":", lw=0.4, alpha=0.6)
        ax.legend(fontsize=9, loc="best")
        ax.invert_xaxis()          # 步长向右增大 = 网格更粗
        fig.tight_layout()
        out = os.path.join(RESULTS, png_name)
        fig.savefig(out, dpi=170)
        plt.close(fig)
        print(f"wrote {out}")


if __name__ == "__main__":
    main()
