"""Exact check of the simulator: expected ON periods from one seed on a short chain, solved over all
2^N - 1 non-dark states, compared with cp_sim (build it first in ../simulation).

Usage:  python3 exact_chain.py
"""
import itertools, subprocess, os
import numpy as np

SIM = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "simulation", "cp_sim")

def exact(L, rate, start):
    states = [s for s in itertools.product((0, 1), repeat=L) if any(s)]
    idx = {s: k for k, s in enumerate(states)}
    A = np.zeros((len(states), len(states))); b = np.zeros(len(states))
    for s in states:
        k = idx[s]; ev = []
        for i in range(L):
            if s[i]:
                ev.append((1.0, i, 0))
            else:
                w = sum(rate[j][i] for j in range(L) if s[j])
                if w > 0: ev.append((w, i, 1))
        tot = sum(e[0] for e in ev); A[k, k] = 1.0
        for w, i, new in ev:
            t = list(s); t[i] = new; t = tuple(t); p = w / tot
            if new == 0: b[k] += p
            if any(t): A[k, idx[t]] -= p
    return np.linalg.solve(A, b)[idx[start]]

def simulate(L, r, cut, runs=400000, seed=4242):
    out = "exact_check"
    subprocess.run([SIM, "lat=0", f"Lx={L}", "periodic=0", f"cut={cut}", f"r1={r}", "mode=1", "h=0", "pb=0",
                    "newlat=0", "tmax=1e9", "tmin=0.1", "K=20", f"runs={runs}", f"seed={seed}", f"out={out}"], check=True)
    d = np.loadtxt(out + ".runs"); n = d[:, 2] * d[:, 0]
    for ext in (".runs", ".ts"): os.remove(out + ext)
    return n.mean(), n.std() / np.sqrt(len(n))

if __name__ == "__main__":
    cases = [(2, 1.957, False), (3, 1.2, False), (5, 0.8, False), (5, 0.8, True), (6, 1.6, False)]
    for L, r, tails in cases:
        if tails:
            rate = [[r / abs(i - j) ** 6 if 0 < abs(i - j) <= 3 else 0.0 for j in range(L)] for i in range(L)]
        else:
            rate = [[r if abs(i - j) == 1 else 0.0 for j in range(L)] for i in range(L)]
        start = tuple(1 if i == L // 2 else 0 for i in range(L))
        ex = exact(L, rate, start)
        m, se = simulate(L, r, 3.05 if tails else 1.05)
        print(f"L={L} r={r}{' tails' if tails else ''}: exact {ex:.4f}   simulated {m:.4f} +- {se:.4f}")
