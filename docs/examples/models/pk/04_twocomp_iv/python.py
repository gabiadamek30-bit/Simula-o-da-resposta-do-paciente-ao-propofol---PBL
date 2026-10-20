#!/usr/bin/env python3
"""
Two-Compartment IV Bolus Model - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("Two-Compartment IV Bolus Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # Model parameters
    CL = 5.0   # Clearance (L/h)
    V1 = 10.0  # Central volume (L)
    Q = 10.0   # Inter-compartmental clearance (L/h)
    V2 = 40.0  # Peripheral volume (L)
    Dose = 100.0  # mg

    # Dense early sampling for distribution phase
    saveat = list(np.arange(0.0, 1.0, 0.05)) + list(np.arange(1.0, 48.5, 0.25))

    # Run simulation using the actual API
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_twocomp_iv_bolus(
        cl=CL,
        v1=V1,
        q=Q,
        v2=V2,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=48.0,
        saveat=saveat
    )

    # Extract results
    times = np.array(result["times"])
    conc = np.array(result["observations"]["conc"])

    # Compute metrics
    cmax = np.max(conc)
    tmax = times[np.argmax(conc)]
    auc = np.trapz(conc, times)

    # Derived parameters
    Vss = V1 + V2
    k10 = CL / V1
    k12 = Q / V1
    k21 = Q / V2

    # Macro rate constants (eigenvalues)
    sum_k = k10 + k12 + k21
    diff = np.sqrt((k10 + k12 + k21)**2 - 4*k10*k21)
    alpha = (sum_k + diff) / 2
    beta = (sum_k - diff) / 2

    t_half_alpha = np.log(2) / alpha
    t_half_beta = np.log(2) / beta

    print("\nResults:")
    print("-" * 50)
    print("Micro-constants:")
    print(f"  k10 = {k10:.4f} 1/h")
    print(f"  k12 = {k12:.4f} 1/h")
    print(f"  k21 = {k21:.4f} 1/h")
    print("\nMacro-constants:")
    print(f"  α = {alpha:.4f} 1/h (t½α = {t_half_alpha:.2f} h)")
    print(f"  β = {beta:.4f} 1/h (t½β = {t_half_beta:.2f} h)")
    print("\nVolumes:")
    print(f"  V1 = {V1} L (central)")
    print(f"  V2 = {V2} L (peripheral)")
    print(f"  Vss = {Vss} L (steady-state)")
    print("\nSimulated PK:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  AUC  = {auc:.2f} mg·h/L")

    # Show bi-exponential profile
    print("\nConcentration-Time Profile (bi-exponential):")
    print("-" * 45)
    print("Time (h)  Conc (mg/L)   Phase")
    print("-" * 45)
    for t in [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0]:
        idx = np.argmin(np.abs(times - t))
        phase = "distribution (α)" if t < 1.0 else "terminal (β)"
        print(f"  {t:5.1f}     {conc[idx]:8.4f}    {phase}")

    return result


if __name__ == "__main__":
    main()
