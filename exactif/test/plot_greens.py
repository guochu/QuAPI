#!/usr/bin/env python3
"""exactif/test/plot_greens.py —— Toulouse（U=0）严格核 1/2/4 阶 Trotter
Green 函数误差标度图（数据：greens_error.dat，由 toulouse_greens.jl 生成）。"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

dat = np.loadtxt("greens_error.dat")
fig, ax = plt.subplots(figsize=(6.2, 4.6))
colors = {1: "#c0392b", 2: "#2471a3", 4: "#1e8449", 6: "#8e44ad"}
labels = {1: "Trotter 1st (Lie)", 2: "Trotter 2nd (Strang)", 4: "Trotter 4th (Suzuki-5)", 6: "Trotter 6th (Suzuki-recursive)"}
for order in (1, 2, 4, 6):
    d = dat[dat[:, 0] == order]
    N, dt, err = d[:, 1], d[:, 2], d[:, 3]
    ax.loglog(dt, err, "o-", color=colors[order], label=labels[order], ms=5, lw=1.4)
    # 斜率标注（取对数拟合，显示在每组右侧）
    if len(dt) >= 2:
        sl = np.polyfit(np.log(dt), np.log(err), 1)[0]
        ax.annotate(f"slope ≈ {abs(sl):.2f}", xy=(dt[-1], err[-1]),
                    xytext=(6, 4), textcoords="offset points",
                    color=colors[order], fontsize=9)
ref = np.array([1, 0.025])
for p, c in ((1, "#7f8c8d"), (2, "#95a5a6"), (4, "#bdc3c7"), (6, "#d2b4de")):
    ax.loglog(ref, [3e-2 * (ref[0] / t) ** p for t in ref], "--", color=c, lw=1,
              label=f"Δt^{p} guide")
ax.set_xlabel(r"$\Delta t$")
ax.set_ylabel(r"$|G^>_{\rm trot}(T) - G^>_{\rm exact}(T)|$")
ax.set_title("Toulouse (U=0) Green function via exact kernel\n"
             r"$T=2.0,\ \beta=2.0,\ \epsilon_d=1.0,\ \Gamma_0=0.15$, GL-64 bath, machine-precision kernel data")
ax.legend(fontsize=8, loc="lower right")
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig("greens_error.png", dpi=160)
print("saved greens_error.png")
