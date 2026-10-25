"""Burst statistics on small lattices below threshold.
(1) 'true' single-seed avalanches (h=0, random seed site)  vs  (2) bursts as an experiment would see them
(dark lattice, spontaneous wake-ups at rate h, bursts = stretches with at least one dye ON)."""
import os; os.makedirs("out", exist_ok=True)
import subprocess, numpy as np, json, os
RS=[0.0,0.1,0.2,0.25,0.3,0.35,0.4,0.5]
res=[]
for L in (6,12,24):
    for pbn in (15,50):
        for r in RS:
            # (1) single-seed, pristine lattice (occupancy 0.85)
            out=f"out/av_seed_{L}_{pbn}_{r}"
            subprocess.run(["./cp_sim","lat=3",f"Lx={L}",f"Ly={L}","periodic=0","cut=3.05",f"r1={r}","mode=1","seedrand=1","h=0",f"pb={1/pbn}",
                            "occ=0.85","tmax=1e6","tmin=0.1","K=50","runs=20000","seed=77",f"out={out}"],check=True)
            d=np.loadtxt(out+".runs"); S=d[:,1]*d[:,0]; n=d[:,2]*d[:,0]
            rec=dict(L=L,pbn=pbn,r1=r,seed_meanS=float(S.mean()),seed_meanN=float(n.mean()),seed_cv=float(S.std()/S.mean()),
                     seed_p99=float(np.percentile(n,99)),seed_max=float(n.max()))
            os.remove(out+".runs"); os.remove(out+".ts")
            for h in (0.002,0.0005):
                out=f"out/av_field_{L}_{pbn}_{r}_{h}"
                subprocess.run(["./cp_sim","lat=3",f"Lx={L}",f"Ly={L}","periodic=0","cut=3.05",f"r1={r}","mode=2",f"h={h}",f"pb={1/pbn}",
                                "occ=0.85","tmax=1e7","tmin=0.1","K=50","runs=400","seed=78","fdead=0.05","maxbursts=4000000",f"out={out}"],check=True)
                b=np.fromfile(out+".bursts").reshape(-1,5)
                rec[f"h{h}_bursts_per_lattice"]=len(b)/400
                rec[f"h{h}_meanS"]=float(b[:,0].mean()); rec[f"h{h}_meanN"]=float(b[:,2].mean()); rec[f"h{h}_meanT"]=float(b[:,1].mean())
                os.remove(out+".bursts"); os.remove(out+".runs"); os.remove(out+".ts")
            res.append(rec); print(rec,flush=True)
os.makedirs("results", exist_ok=True)
json.dump(res,open("results/aval_results.json","w"))
