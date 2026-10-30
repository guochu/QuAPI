#!/usr/bin/env python3
"""虚时版误差标度图（greens_imagtime_error.dat）。"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dat = np.loadtxt("greens_imagtime_error.dat")
fig, ax = plt.subplots(figsize=(6.2, 4.6))
colors = {1: "#c0392b", 2: "#2471a3", 4: "#1e8449", 6: "#8e44ad"}
labels = {1: "Trotter 1st", 2: "Trotter 2nd", 4: "Trotter 4th", 6: "Trotter 6th"}
for order in (1, 2, 4, 6):
    d = dat[dat[:, 0] == order]
    ax.loglog(d[:, 1], d[:, 2], "o-", color=colors[order], label=labels[order], ms=5, lw=1.4)
    if len(d) >= 2:
        sl = np.polyfit(np.log(d[:, 1]), np.log(d[:, 2]), 1)[0]
        ax.annotate(f"slope ≈ {abs(sl):.2f}", xy=(d[-1, 1], d[-1, 2]),
                    xytext=(6, 4), textcoords="offset points",
                    color=colors[order], fontsize=9)
for p, c in ((2, "#95a5a6"), (4, "#bdc3c7"), (6, "#d2b4de")):
    xs = np.array([5, 40.0])
    ax.loglog(xs, 3e-2 * (xs[0] / xs) ** p, "--", color=c, lw=1, label=f"$N_\\tau^{{-{p}}}$ guide")
ax.axhline(3e-13, color="k", ls=":", lw=1)
ax.annotate("double-precision floor", xy=(6, 3e-13), xytext=(0, 8),
            textcoords="offset points", fontsize=8)
ax.set_xlabel(r"$N_\tau$ (at fixed $\beta=2.0$)")
ax.set_ylabel(r"$|G(\beta/2)_{\rm trot} - G(\beta/2)_{\rm exact}|$")
ax.set_title("Toulouse (U=0) imaginary-time Green function via exact kernel\n"
             r"$\beta=2.0,\ \epsilon_d=1.0,\ \Gamma_0=0.15$, GL-64 bath, machine-precision kernel data")
ax.legend(fontsize=8, loc="lower left")
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("greens_imagtime_error.png", dpi=160)
print("saved greens_imagtime_error.png")
