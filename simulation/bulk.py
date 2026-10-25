"""Bulk (periodic, 256x256) runs: (W) critical quench with field / bleaching, (T) thresholds with vacancies
and for the denser every-helix grid.  Results -> results/*.ts (columns: t, ON fraction, <rho^2>, survival, <N_on>, R^2, cumulative ON-time per site, alive fraction)."""
import os; os.makedirs("results", exist_ok=True)
import subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
base = dict(lat=3, Lx=256, Ly=256, periodic=1, cut=3.05, mode=0, tmax=1000, tmin=0.1, K=140, runs=2)
jobs = []
# W: scaling window at the clean critical point of the pixel grid (r1c ~ 0.334)
for tag, extra in [("clean", {}), ("h2e-3", dict(h=0.002)), ("h5e-4", dict(h=0.0005)), ("h1e-4", dict(h=0.0001)),
                   ("pb15", dict(pb=1/15)), ("pb50", dict(pb=1/50)), ("pb145", dict(pb=1/145)),
                   ("h2e-3_pb50", dict(h=0.002, pb=1/50))]:
    jobs.append((f"results/W_{tag}", dict(base, r1=0.334, seed=101, **extra)))
# T1: threshold of the pixel grid with random vacancies (no field, no bleaching)
for occ, rs in [(0.85, [0.38, 0.39, 0.40, 0.42, 0.44, 0.46]), (0.70, [0.46, 0.49, 0.52, 0.56, 0.60, 0.64])]:
    for r in rs:
        jobs.append((f"results/T_occ{occ}_{r}", dict(base, r1=r, occ=occ, seed=202)))
# T2: denser grid, a dye on every helix (rows 2.65 nm apart): P=10.88, b=2.65 -> nn 6.05 nm
for r in [0.14, 0.15, 0.16, 0.17]:
    jobs.append((f"results/T_dense_{r}", dict(base, r1=r, b=2.65, seed=303)))
# T3: denser grid at the fitted coupling r(6.05 nm)=1.91, thinned to occupancy p
for occ in [0.10, 0.14, 0.18, 0.24]:
    jobs.append((f"results/T_densefit_occ{occ}", dict(base, r1=1.91, b=2.65, occ=occ, seed=404)))
# S: single seed at the clean critical point (survival, activity, spreading radius); lattice built once
jobs.append(("results/S_clean", dict(lat=3, Lx=96, Ly=96, periodic=1, cut=3.05, r1=0.334, mode=1, newlat=0, tmax=100, tmin=0.5, K=50, runs=100000, seed=909)))
def run(j):
    out, a = j
    subprocess.run(["./cp_sim"] + [f"{k}={v}" for k, v in a.items()] + [f"out={out}"], check=True)
    return out
t0 = time.time()
with ThreadPoolExecutor(2) as ex:
    for o in ex.map(run, jobs):
        print("done", o, f"{time.time()-t0:.0f}s", flush=True)
print("ALLDONE", flush=True)
