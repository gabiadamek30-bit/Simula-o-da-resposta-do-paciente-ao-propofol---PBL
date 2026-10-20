#!/usr/bin/env python3
"""
One-Compartment IV Bolus Model - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("One-Compartment IV Bolus Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # Model parameters
    CL = 5.0   # Clearance (L/h)
    V = 50.0   # Volume of distribution (L)
    Dose = 100.0  # mg

    # Run simulation using the actual API
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_iv_bolus(
        cl=CL,
        v=V,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=48.0,
        saveat=[float(t) for t in np.arange(0.0, 48.5, 0.5)]
    )

    # Extract results
    times = np.array(result["times"])
    conc = np.array(result["observations"]["conc"])

    # Compute metrics
    cmax = np.max(conc)
    tmax = times[np.argmax(conc)]
    auc = np.trapz(conc, times)

    # Theoretical values
    k = CL / V
    t_half_theoretical = np.log(2) / k
    cmax_theoretical = Dose / V
    auc_theoretical = Dose / CL

    print("\nResults:")
    print("-" * 50)
    print("Simulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Tmax = {tmax:.2f} h")
    print(f"  AUC  = {auc:.2f} mg·h/L")
    print("\nTheoretical:")
    print(f"  Cmax = {cmax_theoretical:.4f} mg/L")
    print(f"  t½   = {t_half_theoretical:.2f} h")
    print(f"  AUC  = {auc_theoretical:.2f} mg·h/L")

    # Sample output
    print("\nConcentration-Time Profile (first 12h):")
    print("-" * 30)
    for t in range(0, 13, 2):
        idx = np.where(np.isclose(times, t))[0]
        if len(idx) > 0:
            print(f"t={t:2d}h: {conc[idx[0]]:.4f} mg/L")

    return result


if __name__ == "__main__":
    main()
