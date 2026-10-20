#!/usr/bin/env python3
"""
One-Compartment Oral First-Order Absorption - Python Example

Run: python python.py
"""

import neopkpd
import numpy as np


def main():
    print("One-Compartment Oral First-Order Absorption Model")
    print("=" * 50)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # Model parameters
    Ka = 1.5   # Absorption rate constant (1/h)
    CL = 5.0   # Clearance (L/h)
    V = 50.0   # Volume of distribution (L)
    Dose = 100.0  # mg

    # Run simulation using the actual API
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_oral_first_order(
        ka=Ka,
        cl=CL,
        v=V,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=48.0,
        saveat=[float(t) for t in np.arange(0.0, 48.25, 0.25)]
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
    tmax_theoretical = np.log(Ka/k) / (Ka - k)
    cmax_theoretical = (Dose/V) * (Ka/(Ka-k)) * (np.exp(-k*tmax_theoretical) - np.exp(-Ka*tmax_theoretical))

    print("\nResults:")
    print("-" * 50)
    print(f"Parameters: Ka={Ka} 1/h, CL={CL} L/h, V={V} L")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Tmax = {tmax:.2f} h")
    print(f"  AUC  = {auc:.2f} mg·h/L")
    print("\nTheoretical:")
    print(f"  Cmax = {cmax_theoretical:.4f} mg/L")
    print(f"  Tmax = {tmax_theoretical:.2f} h")

    # Absorption vs elimination phases
    print("\nConcentration-Time Profile:")
    print("-" * 40)
    print("Time (h)  Conc (mg/L)  Phase")
    print("-" * 40)
    for t in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]:
        idx = np.argmin(np.abs(times - t))
        if t < tmax:
            phase = "absorption"
        elif np.isclose(t, tmax, atol=0.25):
            phase = "peak"
        else:
            phase = "elimination"
        print(f"  {t:4.1f}     {conc[idx]:8.4f}    {phase}")

    return result


if __name__ == "__main__":
    main()
