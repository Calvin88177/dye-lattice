"""Figures for dp_forster_theory.tex. Reads the critical-window runs in ../simulation/results."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter, LogLocator

RES = "../simulation/results/"
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE3 = ["#86b6ef", "#256abf", "#0d366b"]
ACC, ORANGE, GREY = "#2a78d6", "#eb6834", "#898781"
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.03})
DELTA, NPAR, SIG = 0.4505, 1.295, 2.18
EH, EB = NPAR / SIG, 1 / (1 - DELTA + 1 / NPAR)

def ts(tag):
    d = np.loadtxt(RES + f"W_{tag}.ts"); return d[:, 0], d[:, 1]

# ---- Fig 1: coupling against distance -------------------------------------------------
d = np.linspace(2, 12, 400)
B = lambda d, d0: 1 / (1 + (d / d0) ** 6)
fig, ax = plt.subplots(figsize=(5.4, 3.0))
ax.plot(d, B(d, 6.74), color=ACC, lw=1.8, label=r"fitted, $d_0=6.74$ nm")
ax.plot(d, B(d, 2.38), color=GREY, lw=1.8, label=r"weak end, $d_0=2.38$ nm")
for bc, lab in [(0.61, "chain"), (0.28, "pixel grid, 85% labelled"), (0.15, "every-helix grid")]:
    ax.axhline(bc, color=INK2, lw=0.7, ls=(0, (4, 3)))
    ax.text(11.9, bc + 0.015, f"{lab}: $B_c={bc}$", ha="right", va="bottom", fontsize=7.5, color=INK2)
for dd, lab in [(6.05, "6.05"), (7.59, "7.59")]:
    ax.plot([dd], [B(dd, 6.74)], "o", ms=5, color=ACC, mec="white", mew=1)
    ax.annotate(f"{lab} nm", (dd, B(dd, 6.74)), xytext=(-7, 5), textcoords="offset points", ha="right", va="bottom", fontsize=7.5, color=INK)
ax.set_xlim(2, 12); ax.set_ylim(0, 1.02)
ax.set_xlabel("dye spacing $d$ (nm)"); ax.set_ylabel("hand-over probability $B$")
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 0.98))
fig.savefig("fig_coupling.pdf"); plt.close(fig)

# ---- Fig 2: critical decay and its cut-offs ---------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.7), sharey=True)
tc, rc = ts("clean")
for ax, items, title, expo in [
        (axs[0], [("h1e-4", 1e-4), ("h5e-4", 5e-4), ("h2e-3", 2e-3)], "spontaneous recovery $h$", EH),
        (axs[1], [("pb145", 1 / 145), ("pb50", 1 / 50), ("pb15", 1 / 15)], r"bleaching $p_\mathrm{b}$", EB)]:
    m = tc <= 1000
    ax.plot(tc[m], rc[m] * tc[m] ** DELTA, color=GREY, lw=1.4, label="unperturbed")
    for (tag, s), col in zip(items, BLUE3):
        t, r = ts(tag); m = (t <= 1000) & (r > 0)
        lab = (f"$h={s:g}$" if tag.startswith("h") else rf"$1/p_\mathrm{{b}}={round(1/s)}$")
        ax.plot(t[m], r[m] * t[m] ** DELTA, color=col, lw=1.4, label=lab)
        tp = s ** -expo
        ax.plot([tp], [0.04], marker="v", ms=6, color=col, mec="white", mew=0.8, clip_on=False)
    ax.set_xscale("log"); ax.set_xlim(0.3, 1000); ax.set_ylim(0, 1.6)
    ax.set_title(title, fontsize=9, color=INK, loc="left")
    ax.set_xlabel("time $t$ (ON lifetimes)")
    ax.legend(loc="upper left", fontsize=7.2)
axs[0].set_ylabel(r"$\rho(t)\,t^{\delta}$")
axs[1].text(0.98, 0.92, "▼ predicted crossover", transform=axs[1].transAxes, ha="right", fontsize=7.2, color=INK2)
fig.savefig("fig_window.pdf"); plt.close(fig)

# ---- Fig 3: crossover times against strength -------------------------------------------
def crossover(tag, frac):
    t, r = ts(tag); r0 = np.interp(t, tc, rc)
    dev = np.where(np.abs(r / r0 - 1) > frac)[0]; return t[dev[0]]
fig, axs = plt.subplots(1, 2, figsize=(6.6, 2.6))
hs = np.array([1e-4, 5e-4, 2e-3]); hsim = [crossover(t, 0.2) for t in ("h1e-4", "h5e-4", "h2e-3")]
pbs = np.array([1 / 145, 1 / 50, 1 / 15]); bsim = [crossover(t, 0.5) for t in ("pb145", "pb50", "pb15")]
for ax, x, y, e, xl, note in [(axs[0], hs, hsim, EH, "$h$ (per ON lifetime)", r"$t_h\sim h^{-%.2f}$" % EH),
                              (axs[1], pbs, bsim, EB, r"$p_\mathrm{b}$ (per ON period)", r"$t_b\sim p_\mathrm{b}^{-%.2f}$" % EB)]:
    xx = np.logspace(np.log10(x.min() / 2), np.log10(x.max() * 2), 50)
    ax.plot(xx, xx ** -e, color=ACC, lw=1.6, label="prediction (prefactor 1)")
    ax.plot(x, y, "o", ms=6, color=INK, mec="white", mew=1, label="simulation")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.yaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(xl); ax.text(0.95, 0.88, note, transform=ax.transAxes, ha="right", fontsize=9, color=INK)
axs[0].set_ylabel("crossover time (ON lifetimes)")
axs[0].legend(loc="lower left", fontsize=7.2, bbox_to_anchor=(0.0, 0.0))
axs[0].set_title("20% departure", fontsize=8, color=INK2, loc="left"); axs[1].set_title("50% departure", fontsize=8, color=INK2, loc="left")
fig.savefig("fig_crossover.pdf"); plt.close(fig)
print("h sim", hsim, "pb sim", bsim)
