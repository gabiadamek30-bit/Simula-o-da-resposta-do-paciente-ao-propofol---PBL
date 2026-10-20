#!/usr/bin/env python3
"""
Two-Compartment Oral Model - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("Two-Compartment Oral Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # Model parameters
    Ka = 1.5   # Absorption rate constant (1/h)
    CL = 5.0   # Clearance (L/h)
    V1 = 10.0  # Central volume (L)
    Q = 10.0   # Inter-compartmental clearance (L/h)
    V2 = 40.0  # Peripheral volume (L)
    Dose = 100.0  # mg

    # Run simulation using the actual API
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_twocomp_oral(
        ka=Ka,
        cl=CL,
        v1=V1,
        q=Q,
        v2=V2,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=72.0,
        saveat=[float(t) for t in np.arange(0.0, 72.25, 0.25)]
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

    print("\nResults:")
    print("-" * 50)
    print("Parameters:")
    print(f"  Ka = {Ka} 1/h, CL = {CL} L/h")
    print(f"  V1 = {V1} L, V2 = {V2} L, Q = {Q} L/h")
    print(f"  Vss = {Vss} L")
    print("\nSimulated PK:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Tmax = {tmax:.2f} h")
    print(f"  AUC  = {auc:.2f} mg·h/L")

    # Show triphasic profile
    print("\nConcentration-Time Profile:")
    print("-" * 45)
    for t in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0, 72.0]:
        idx = np.argmin(np.abs(times - t))
        print(f"t={t:5.1f}h: {conc[idx]:.4f} mg/L")

    return result


if __name__ == "__main__":
    main()
