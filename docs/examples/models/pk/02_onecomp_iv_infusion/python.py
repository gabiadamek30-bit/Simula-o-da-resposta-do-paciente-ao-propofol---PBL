#!/usr/bin/env python3
"""
One-Compartment IV Infusion Model - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("One-Compartment IV Infusion Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # Model parameters
    CL = 5.0        # Clearance (L/h)
    V = 50.0        # Volume of distribution (L)
    Dose = 100.0    # mg
    Duration = 1.0  # Infusion duration (h)

    # Dense sampling during and after infusion
    saveat = list(np.arange(0.0, 2.1, 0.1)) + list(np.arange(2.5, 24.5, 0.5))

    # Run simulation with infusion (use duration parameter)
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_iv_bolus(
        cl=CL,
        v=V,
        doses=[{"time": 0.0, "amount": Dose, "duration": Duration}],
        t0=0.0,
        t1=24.0,
        saveat=saveat
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
    Rate = Dose / Duration
    C_end_infusion = (Rate / CL) * (1 - np.exp(-k * Duration))

    print("\nResults:")
    print("-" * 50)
    print(f"Infusion: {Dose} mg over {Duration} h (Rate = {Rate} mg/h)")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Tmax = {tmax:.2f} h (end of infusion)")
    print(f"  AUC  = {auc:.2f} mg·h/L")
    print("\nTheoretical C at end of infusion:")
    print(f"  C(Tinf) = {C_end_infusion:.4f} mg/L")

    # Compare to bolus
    C_bolus_max = Dose / V
    print("\nComparison to IV bolus:")
    print(f"  Bolus Cmax = {C_bolus_max:.4f} mg/L")
    print(f"  Infusion Cmax/Bolus Cmax = {cmax/C_bolus_max * 100:.1f}%")

    # Sample output
    print("\nConcentration during and after infusion:")
    print("-" * 35)
    for t in [0.0, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0, 12.0]:
        idx = np.argmin(np.abs(times - t))
        phase = "(infusion)" if t <= Duration else "(post-infusion)"
        print(f"t={t:4.1f}h: {conc[idx]:.4f} mg/L {phase}")

    return result


if __name__ == "__main__":
    main()
