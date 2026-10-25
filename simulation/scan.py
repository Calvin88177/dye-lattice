"""Threshold scans: prints rho(t)*t^delta at several times; flat = critical, rising = active, falling = absorbing.
Usage: python3 scan.py validate | chain_tail | tri_tail | pix_tail"""
import subprocess, sys, numpy as np, itertools, os
os.makedirs("out", exist_ok=True)
from concurrent.futures import ThreadPoolExecutor
def run(args, out):
    cmd=["./cp_sim"]+[f"{k}={v}" for k,v in args.items()]+[f"out={out}"]
    subprocess.run(cmd, check=True)
    return out
def scan(tag, base, rs, delta, times, workers=2):
    jobs=[]
    with ThreadPoolExecutor(workers) as ex:
        futs=[ex.submit(run, dict(base, r1=r), f"out/{tag}_{r}") for r in rs]
        outs=[f.result() for f in futs]
    for r,o in zip(rs,outs):
        d=np.loadtxt(o+".ts"); t,rho=d[:,0],d[:,1]
        sel=[np.argmin(abs(t-x)) for x in times]
        print(tag, r, " ".join(f"{rho[i]*t[i]**delta:.4f}" for i in sel), flush=True)
if __name__=="__main__":
    which=sys.argv[1]
    if which=="chain_tail":
        scan("c1dT", dict(lat=0,Lx=20000,periodic=1,cut=3.05,mode=0,tmax=3000,tmin=0.1,K=120,runs=3,seed=21), [1.56,1.58,1.60,1.62], 0.1595, (10,30,100,300,1000,3000))
    if which=="tri_tail":
        scan("triT", dict(lat=2,Lx=256,Ly=256,periodic=1,cut=3.05,mode=0,tmax=1000,tmin=0.1,K=120,runs=2,seed=22), [0.229,0.232,0.234,0.236,0.242], 0.4505, (10,30,100,300,1000))
    if which=="pix_tail":
        scan("pixT", dict(lat=3,Lx=256,Ly=256,periodic=1,cut=3.05,mode=0,tmax=1000,tmin=0.1,K=120,runs=2,seed=23), [0.33,0.334,0.338,0.342,0.35], 0.4505, (10,30,100,300,1000))
    if which=="validate":   # known critical points: chain 1.6489, square 0.41219, triangular 0.25797 (nearest neighbours only)
        scan("v1d", dict(lat=0,Lx=20000,periodic=1,cut=1.05,mode=0,tmax=3000,tmin=0.1,K=120,runs=3,seed=7), [1.62,1.6489,1.68], 0.1595, (10,30,100,300,1000,3000))
        scan("v2sq", dict(lat=1,Lx=400,Ly=400,periodic=1,cut=1.05,mode=0,tmax=1500,tmin=0.1,K=120,runs=2,seed=11), [0.405,0.41219,0.420], 0.4505, (10,30,100,300,1000))
        scan("v2tri", dict(lat=2,Lx=400,Ly=400,periodic=1,cut=1.05,mode=0,tmax=1500,tmin=0.1,K=120,runs=2,seed=13), [0.2535,0.257967,0.2629], 0.4505, (10,30,100,300,1000))
