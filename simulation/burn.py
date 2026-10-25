"""Burn-curve scan: finite pixel-grid lattices, all dyes ON at t=0, realistic imperfections.
Writes results/burn_results_compact.json: for each scenario, size and coupling r1 -> statistics over lattices.
flash_mean = share of all light emitted before tstar = 4 x (ON periods per dye) lifetimes;
survivors_mean = unbleached dyes at tstar / lattice sites.
Time unit = one mean ON period (1/k_off).  r1 = wake rate per ON nearest neighbour / k_off."""
import os; os.makedirs("out", exist_ok=True)
import subprocess, json, numpy as np, sys, os, time
from concurrent.futures import ThreadPoolExecutor
RS = [0.0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.85, 1.0, 1.3]
SIZES = [6, 12, 24, 48]
COND = {
  "ideal":       dict(h=0.002,  pb=1/15,  occ=1.0,  sigma=0.0),
  "occ85":       dict(h=0.002,  pb=1/15,  occ=0.85, sigma=0.0),
  "occ70":       dict(h=0.002,  pb=1/15,  occ=0.70, sigma=0.0),
  "occ85_jit":   dict(h=0.002,  pb=1/15,  occ=0.85, sigma=0.13),
  "occ85_lowh":  dict(h=0.0005, pb=1/15,  occ=0.85, sigma=0.0),
  "ideal_pb50":  dict(h=0.002,  pb=1/50,  occ=1.0,  sigma=0.0),
  "occ85_pb50":  dict(h=0.002,  pb=1/50,  occ=0.85, sigma=0.0),
  "occ70_pb50":  dict(h=0.002,  pb=1/50,  occ=0.70, sigma=0.0),
  "occ85_pb145": dict(h=0.002,  pb=1/145, occ=0.85, sigma=0.0),
}
def job(cn, L, r):
    c = COND[cn]
    out = f"out/burn_{cn}_{L}_{r}"
    runs = int((400 if L <= 24 else 200) * min(1.0, 15 * c["pb"] * 2))
    runs = max(runs, 40)
    tmax = 40.0 / (c["pb"] * c["h"])          # long enough for every dye to finish
    tstar = 4.0 / c["pb"]
    args = dict(lat=3, Lx=L, Ly=L, periodic=0, cut=3.05, r1=r, mode=0, tmax=tmax, tmin=0.02, K=400,
                runs=runs, seed=1000 + 17 * L + int(1000 * r), newlat=1, tstar=tstar, **c)
    subprocess.run(["./cp_sim"] + [f"{k}={v}" for k, v in args.items()] + [f"out={out}"], check=True)
    d = np.loadtxt(out + ".runs")
    ts = np.loadtxt(out + ".ts")
    os.remove(out + ".runs"); os.remove(out + ".ts")
    f = d[:,10]; a = d[:,11]*d[:,0]/(L*L); r4 = lambda v: [float(f"{x:.4g}") for x in v[::2]]
    return dict(cond=cn, L=L, r1=r, lattices=runs, dyes_present=round(float(d[:,0].mean()),2), on_periods_per_dye=round(float(d[:,2].mean()),3),
                tstar=tstar, flash_mean=round(float(f.mean()),5), flash_sd=round(float(f.std()),5),
                survivors_mean=round(float(a.mean()),5), survivors_sd=round(float(a.std()),5),
                half_burn_time_median=float(f"{np.median(d[:,4]):.5g}"),
                t=r4(ts[:,0]), rho=r4(ts[:,1]), cum=r4(ts[:,6]), alive=r4(ts[:,7]))
if __name__ == "__main__":
    conds = sys.argv[1:] or list(COND)
    t0 = time.time()
    jobs = [(cn, L, r) for cn in conds for L in SIZES for r in RS]
    with ThreadPoolExecutor(2) as ex:
        res = list(ex.map(lambda a: job(*a), jobs))
    os.makedirs("results", exist_ok=True)
    fn = "results/burn_results_compact.json"
    old = json.load(open(fn)) if os.path.exists(fn) else []
    old = [o for o in old if o["cond"] not in conds] + res
    json.dump(old, open(fn, "w"), separators=(",", ":"))
    print("done", conds, f"{time.time()-t0:.0f}s")
