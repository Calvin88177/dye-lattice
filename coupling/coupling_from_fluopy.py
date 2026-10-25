"""Recompute the pair hand-over chance B(d) from the Cy5 spectra and rate constants of the
published model (Ebert, Sauer and Doose, PLOS Comput. Biol. 2026; code at github.com/super-resolution/fluopy).

Usage:
    git clone https://github.com/super-resolution/fluopy
    python3 coupling_from_fluopy.py fluopy/src/fluopy/fluorophore_spectra/cy5_data

Result: R0 = 4.29 nm for transfer from an ON dye to a dark one, d0 = R0 (eta/p0)^(1/6) = 6.77 nm,
and B(d) = 1/(1 + (d/d0)^6), matching the full model (d0 = 6.74 nm) to within 0.001 in B.
"""
import csv, sys
import numpy as np

def load(path):
    rows = []
    with open(path) as f:
        for row in csv.reader(f):
            try:
                rows.append([float(x) for x in row[:2]])
            except ValueError:
                pass
    a = np.array(rows)
    return a[:, 0], a[:, 1]

def main(spec_dir):
    wl_e, em = load(f"{spec_dir}/emission.csv")
    wl_o, eps_off = load(f"{spec_dir}/absorption_off.csv")
    wl = np.arange(max(wl_e.min(), wl_o.min()), min(wl_e.max(), wl_o.max()) + 1e-9, 1.0)
    E = np.interp(wl, wl_e, em); A = np.interp(wl, wl_o, eps_off)
    J = np.trapezoid(E / np.trapezoid(E, wl) * A * wl**4, wl)          # M^-1 cm^-1 nm^4

    # Cy5 on DNA, as in fluopy (cy5_dna) and the paper's dSTORM parameter set
    QY, tau, n, k2 = 0.27, 1.7e-9, 1.33, 2 / 3
    k_isc, k_ts, k_iso = 8.3e5, 5e3, 4e6
    kS_mol, kT_mol, phi = 1e9, 1e8, 1e-3           # thiolate electron transfer, success yield
    eta = 1e-4                                      # fitted wake-up yield per transfer
    conc_mM, pH, pKa = 100.0, 7.5, 9.0              # MEA

    R0 = 0.02108 * (k2 * QY * J / n**4) ** (1 / 6)   # nm
    c = 10 ** (pH - pKa) / (10 ** (pH - pKa) + 1) * conc_mM * 1e-3   # thiolate, M
    kD = 1 / tau
    k_tot = kD + kS_mol * c + k_isc + k_iso
    p0 = phi * kS_mol * c / k_tot + (k_isc / k_tot) * phi * kT_mol * c / (kT_mol * c + k_ts)
    d0 = R0 * (eta / p0) ** (1 / 6)

    print(f"J = {J:.3e} M^-1 cm^-1 nm^4,  R0 = {R0:.2f} nm")
    print(f"thiolate = {c*1e3:.2f} mM,  p0 = {p0:.2e} per excitation,  d0 = {d0:.2f} nm")
    print(f"weak end (eta = 2e-7): d0 = {R0 * (1.95e-7 / p0) ** (1/6):.2f} nm")
    print("\n d (nm)   B fitted   B weak")
    for d in (2.95, 5.20, 6.03, 7.14, 9.51, 17.78):
        Bf = 1 / (1 + (d / d0) ** 6); Bw = 1 / (1 + (d / (R0 * (1.95e-7 / p0) ** (1/6))) ** 6)
        print(f" {d:6.2f}   {Bf:8.3f}   {Bw:8.4f}")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "fluopy/src/fluopy/fluorophore_spectra/cy5_data")
