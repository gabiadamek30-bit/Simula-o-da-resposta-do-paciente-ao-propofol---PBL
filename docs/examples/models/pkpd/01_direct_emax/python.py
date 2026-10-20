#!/usr/bin/env python3
"""
Direct Emax PD Model - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("Direct Emax PD Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # PK parameters
    CL = 5.0   # L/h
    V = 50.0   # L
    Dose = 100.0  # mg

    # PD parameters
    E0 = 0.0      # Baseline effect
    Emax = 100.0  # Maximum effect
    EC50 = 1.0    # mg/L

    # Run PKPD simulation using the actual API
    print("\nRunning PKPD simulation...")
    result = neopkpd.simulate_pkpd_direct_emax(
        cl=CL,
        v=V,
        e0=E0,
        emax=Emax,
        ec50=EC50,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=48.0,
        saveat=[float(t) for t in np.arange(0.0, 48.25, 0.25)]
    )

    # Extract results
    times = np.array(result["times"])
    conc = np.array(result["observations"]["conc"])
    effect = np.array(result["observations"]["effect"])

    # Compute metrics
    cmax = np.max(conc)
    emax_observed = np.max(effect)
    e_at_ec50 = E0 + Emax / 2  # Theoretical effect at EC50

    print("\nResults:")
    print("-" * 50)
    print(f"PK Parameters: CL={CL} L/h, V={V} L")
    print(f"PD Parameters: E0={E0}, Emax={Emax}, EC50={EC50} mg/L")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Max Effect = {emax_observed:.2f}")
    print(f"  Effect at C=EC50 (theoretical) = {e_at_ec50:.2f}")

    # Show PK/PD relationship
    print("\nConcentration-Effect Profile:")
    print("-" * 45)
    print("Time (h)  Conc (mg/L)  Effect")
    print("-" * 45)
    for t in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0]:
        idx = np.argmin(np.abs(times - t))
        c = conc[idx]
        e = effect[idx]
        print(f"  {t:4.0f}     {c:7.3f}    {e:.2f}")

    # Demonstrate no hysteresis
    print("\nNote: Direct effect model - no hysteresis between PK and PD")
    print("Effect at any time = E0 + Emax × C / (EC50 + C)")

    return result


if __name__ == "__main__":
    main()
