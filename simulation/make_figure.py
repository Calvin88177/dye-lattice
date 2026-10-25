"""Six-panel summary figure for the dye-lattice contact-process simulations.
Run from simulation/. Reads results/burn_results_compact.json, results/aval_results.json and results/W_*.ts (written by burn.py, aval.py, bulk.py)."""
import json, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

# ---- palette (validated: ordinal blue ramps, categorical slots 1-3) ----
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE5 = ["#86b6ef", "#5598e7", "#256abf", "#184f95", "#0d366b"]
BLUE4 = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]
BLUE3 = ["#86b6ef", "#256abf", "#0d366b"]
CAT = ["#2a78d6", "#eb6834", "#1baf7a"]
GRAY = "#898781"
LW = 1.5

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.facecolor": SURF, "figure.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.minor.size": 0, "ytick.minor.size": 0,
})

R = json.load(open("results/burn_results_compact.json"))
D = {(o["cond"], o["L"], o["r1"]): o for o in R}
RS = sorted({o["r1"] for o in R})
AV = json.load(open("results/aval_results.json"))

import textwrap
fig, axs = plt.subplots(2, 3, figsize=(15.5, 10.2))
fig.subplots_adjust(left=0.055, right=0.985, top=0.895, bottom=0.065, wspace=0.27, hspace=0.56)

def title(ax, letter, head, sub):
    """Bold one-line heading plus a wrapped three-line note, both kept inside the panel width."""
    lines = textwrap.wrap(sub, width=66)
    assert len(lines) <= 3, (letter, len(lines))
    lines += [""] * (3 - len(lines))
    ax.text(0.0, 1.225, f"{letter}   {head}", transform=ax.transAxes, fontsize=10.5, fontweight="bold", color=INK, va="bottom", ha="left")
    ax.text(0.0, 1.035, "\n".join(lines), transform=ax.transAxes, fontsize=8.1, color=INK2, va="bottom", ha="left", linespacing=1.3)

# ---------------- (a) scale sketch ----------------
ax = axs[0, 0]
ax.set_aspect("equal"); ax.grid(False)
tile = 88.0
sizes = [(8, "8 × 8 tiles · 8,704 pixels · yield ≈ 4 %"), (4, "4 × 4 · 2,176 · ≈ 46 %"), (2, "2 × 2 · 544 · ≈ 86 %"), (1, "1 tile · 136 pixels")]
for (m, lab), col in zip(sizes, BLUE4):
    s = m * tile
    ax.add_patch(Rectangle((0, 0), s, s, facecolor=col, alpha=0.10, edgecolor="none", zorder=1))
    ax.add_patch(Rectangle((0, 0), s, s, facecolor="none", edgecolor=col, lw=1.4, zorder=2))
ax.text(8 * tile - 8, 8 * tile - 14, sizes[0][1], ha="right", va="top", fontsize=8.3, color=INK)
ax.text(4 * tile + 8, 4 * tile - 6, sizes[1][1], ha="left", va="top", fontsize=8.3, color=INK)
ax.text(2 * tile + 8, 2 * tile - 6, sizes[2][1], ha="left", va="top", fontsize=8.3, color=INK)
ax.text(1 * tile + 8, 0.5 * tile, sizes[3][1], ha="left", va="center", fontsize=8.3, color=INK)
# diffraction-limited spot (FWHM ~ 250 nm)
cx, cy, rad = 505, 480, 125
ax.add_patch(Circle((cx, cy), rad, facecolor=GRAY, alpha=0.22, edgecolor="none", zorder=3))
ax.add_patch(Circle((cx, cy), rad, facecolor="none", edgecolor=INK2, lw=1.0, zorder=4))
ax.text(cx, cy + 32, "one diffraction-\nlimited spot", ha="center", va="center", fontsize=8.3, color=INK, zorder=5)
ax.text(cx, cy - 50, "≈ 250 nm across\n≈ 850 pixels", ha="center", va="center", fontsize=7.8, color=INK2, zorder=5, linespacing=1.25)
# scale bar
ax.plot([40, 240], [-52, -52], color=INK, lw=2, solid_capstyle="butt", clip_on=False)
ax.text(140, -68, "200 nm", ha="center", va="top", fontsize=8.3, color=INK, clip_on=False)
ax.set_xlim(-10, 8 * tile + 10); ax.set_ylim(-110, 8 * tile + 10)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values(): sp.set_visible(False)
title(ax, "a", "The lattice is far below optical resolution",
      "Published tile arrays to scale; pixel pitch about 7.6 nm. At the critical point, activity from one "
      "wake-up spreads only 20–75 nm before bleaching or spontaneous wake-ups end the critical regime.")

# ---------------- (b) the flash ----------------
ax = axs[0, 1]
cond, L = "occ85_pb50", 48
o0 = D[(cond, L, 0.0)]
ax.plot(o0["t"], o0["rho"], color=GRAY, lw=LW, label="0 (uncoupled)")
for r, col in zip([0.2, 0.35, 0.5, 0.7, 1.0], BLUE5):
    o = D[(cond, L, r)]
    ax.plot(o["t"], o["rho"], color=col, lw=LW, label=f"{r:g}")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 3000); ax.set_ylim(5e-4, 1.2)
ax.set_xlabel("time since laser on (ON lifetimes; one ≈ 40 ms at 2.5 kW/cm²)")
ax.set_ylabel("fraction of dyes ON")
leg = ax.legend(title="coupling r₁", loc="lower left", ncol=2, handlelength=1.6, columnspacing=1.2, labelcolor=INK2)
leg.get_title().set_color(INK2); leg.get_title().set_fontsize(8)
ax.annotate("flash: activity sustains\nitself until bleaching\nthins the lattice", xy=(40, 0.33), xytext=(120, 0.55), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), ha="left", va="center")
ax.annotate("slow phase: separate bursts\nfrom spontaneous wake-ups", xy=(700, 0.0085), xytext=(110, 0.075), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), ha="left", va="center")
title(ax, "b", "A movie would show a flash, then slow blinking",
      "Simulated 2,304-site lattice, every dye ON when the laser turns on. 15 % of dyes missing, "
      "50 ON periods per dye before bleaching, spontaneous wake-up 0.002 per lifetime.")

# ---------------- (c) flash fraction vs coupling ----------------
ax = axs[0, 2]
for Ls, col in zip([6, 12, 24, 48], BLUE4):
    y = [D[(cond, Ls, r)]["flash_mean"] for r in RS]
    ax.plot(RS, y, color=col, lw=LW, marker="o", ms=5, mec=SURF, mew=1.0, label=f"{Ls*Ls:,} sites")
ax.set_xlim(-0.03, 1.36); ax.set_ylim(0, 0.66)
ax.axvline(0.40, color=MUTED, lw=0.8, zorder=0)
ax.text(0.415, 0.022, "threshold if dyes never bleached", fontsize=7.8, color=INK2, ha="left", va="bottom")
for x, lab in [(0.001, None), (0.194, "pH 8.0"), (0.49, "pH 7.5"), (1.066, "pH 7.0")]:
    ax.plot([x], [0.648], marker="v", ms=7, color=INK, mec=SURF, mew=0.8, clip_on=False, zorder=6)
    if lab: ax.text(x, 0.622, lab, fontsize=7.8, color=INK, ha="center", va="top")
ax.text(-0.012, 0.590, "weak end", fontsize=7.8, color=INK, ha="left", va="top")
ax.set_xlabel("coupling r₁ (wake rate per ON neighbour ÷ off rate)")
ax.set_ylabel("share of all light emitted in the flash")
ax.legend(title="lattice size", loc="lower right", bbox_to_anchor=(1.0, 0.05), handlelength=1.6, labelcolor=INK2).get_title().set_color(INK2)
title(ax, "c", "The threshold appears as a size-dependent rise",
      "Share of a lattice's total light emitted in the first 200 lifetimes (about 8 s), imperfections as in b. "
      "▼ model values for the 7.6 nm grid: fitted coupling at three pH, and the weak end.")

# ---------------- (d) survivors ----------------
ax = axs[1, 0]
for (c, lab), col in zip([("ideal_pb50", "100 % of sites labelled"), ("occ85_pb50", "85 %"), ("occ70_pb50", "70 %")], CAT):
    y = [D[(c, 48, r)]["survivors_mean"] for r in RS]
    ax.plot(RS, y, color=col, lw=LW, marker="o", ms=5, mec=SURF, mew=1.0, label=lab)
ax.set_xlim(-0.03, 1.36); ax.set_ylim(0.25, 1.05)
ax.set_xlabel("coupling r₁")
ax.set_ylabel("fraction of lattice sites still unbleached after the flash")
ax.legend(title="starting occupancy", loc="upper right", handlelength=1.6, labelcolor=INK2).get_title().set_color(INK2)
ax.annotate("above threshold the lattice burns down\nto the same density wherever it started", xy=(1.0, 0.425), xytext=(0.62, 0.70), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), ha="left", va="center")
title(ax, "d", "Bleaching tunes the lattice to its threshold",
      "Dyes left after the flash (2,304 sites, 50 ON periods per dye) for three starting label densities. "
      "Rydberg gases show the same loss-driven self-organisation (Helmrich et al. 2020).")

# ---------------- (e) scaling window ----------------
ax = axs[1, 1]
def W(tag):
    d = np.loadtxt(f"results/W_{tag}.ts"); return d[:, 0], d[:, 1]
t, r = W("clean"); m = t <= 500
ax.plot(t[m], r[m], color=GRAY, lw=LW, label="ideal: no bleaching, no wake-up")
for tag, lab, col in [("pb145", "145 ON periods per dye", BLUE3[0]), ("pb50", "50", BLUE3[1]), ("pb15", "15", BLUE3[2])]:
    t, r = W(tag); m = (t <= 500) & (r > 0)
    ax.plot(t[m], r[m], color=col, lw=LW, label=lab)
t, r = W("h2e-3"); m = t <= 500
ax.plot(t[m], r[m], color=CAT[1], lw=LW, label="wake-up 0.002 per lifetime, no bleaching")
tg = np.array([3.0, 60.0]); ax.plot(tg, 1.05 * tg ** -0.4505, color=INK, lw=0.8)
ax.text(16, 0.36, "directed percolation, t$^{-0.45}$", fontsize=7.8, color=INK, ha="left", va="bottom", rotation=-21, rotation_mode="anchor")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.3, 500); ax.set_ylim(2e-3, 1.2)
ax.set_xlabel("time since laser on (ON lifetimes)")
ax.set_ylabel("fraction of dyes ON")
ax.legend(loc="lower left", handlelength=1.6, labelcolor=INK2)
title(ax, "e", "The critical power law lasts about one decade",
      "Large lattice (65,536 sites) set at its threshold. Bleaching pulls it off the ideal curve "
      "after 4–20 lifetimes; spontaneous wake-up does so after about 35.")

# ---------------- (f) bursts ----------------
ax = axs[1, 2]
for Ls, col in zip([6, 12, 24], BLUE3):
    rows = [a for a in AV if a["L"] == Ls and a["pbn"] == 50]
    x = [a["r1"] for a in rows]
    ax.plot(x, [a["seed_meanN"] for a in rows], color=col, lw=LW, label=f"{int(round(Ls*Ls*0.85)):,} dyes")
    ax.plot(x, [a["h0.0005_meanN"] for a in rows], ls="none", marker="o", ms=5.5, color=col, mec=SURF, mew=1.0)
ax.set_yscale("log"); ax.set_xlim(-0.02, 0.52); ax.set_ylim(0.8, 3000)
ax.axvline(0.40, color=MUTED, lw=0.8, zorder=0)
ax.text(0.392, 1900, "threshold", fontsize=7.8, color=INK2, ha="right", va="top")
ax.set_xlabel("coupling r₁")
ax.set_ylabel("mean ON periods per burst")
ax.legend(title="lattice size", loc="upper left", handlelength=1.6, labelcolor=INK2).get_title().set_color(INK2)
ax.annotate("dots above the line: bursts\nfrom different seeds merge", xy=(0.35, 27), xytext=(0.03, 110), fontsize=8, color=INK2,
            arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6), ha="left", va="center")
title(ax, "f", "Below threshold, burst size is the observable",
      "Lines: true avalanche from one spontaneous wake-up. Dots: bursts as a movie would record them "
      "(wake-up 0.0005 per lifetime). They agree up to about 120 dyes per lattice.")

fig.savefig("dye_lattice_simulation.png", dpi=150)
print("saved")
