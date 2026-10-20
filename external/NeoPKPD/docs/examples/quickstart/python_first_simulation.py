#!/usr/bin/env python3
"""
NeoPKPD Quickstart - Python

Run:
    source packages/python/.venv/bin/activate
    python docs/examples/quickstart/python_first_simulation.py
"""

import neopkpd
import numpy as np


def main():
    print("NeoPKPD Quickstart - Python")
    print("=" * 40)

    # Initialize Julia backend (required once per session)
    neopkpd.init_julia()

    # 1. Define model parameters
    CL = 5.0   # Clearance (L/h)
    V = 50.0   # Volume of distribution (L)
    Dose = 100.0  # mg

    print(f"\nModel: OneCompIVBolus")
    print(f"Parameters: CL = {CL} L/h, V = {V} L")
    print(f"Dose: {Dose} mg IV at t=0")

    # 2. Run simulation (0-24h with 0.5h intervals)
    print("\nRunning simulation...")
    result = neopkpd.simulate_pk_iv_bolus(
        cl=CL,
        v=V,
        doses=[{"time": 0.0, "amount": Dose}],
        t0=0.0,
        t1=24.0,
        saveat=[float(t) for t in np.arange(0.0, 24.5, 0.5)]
    )

    # 3. Extract results
    times = np.array(result["times"])
    conc = np.array(result["observations"]["conc"])

    # 4. Compute metrics
    cmax = np.max(conc)
    tmax_idx = np.argmax(conc)
    tmax = times[tmax_idx]

    # AUC by trapezoidal rule
    auc = np.trapz(conc, times)

    # Terminal half-life (from last few points)
    log_conc = np.log(conc[-5:])
    t_segment = times[-5:]
    slope = (log_conc[-1] - log_conc[0]) / (t_segment[-1] - t_segment[0])
    t_half = -np.log(2) / slope

    # 5. Display results
    print("\n" + "=" * 40)
    print("RESULTS")
    print("=" * 40)
    print(f"Cmax:     {cmax:.3f} mg/L")
    print(f"Tmax:     {tmax:.1f} h")
    print(f"AUC_0_24: {auc:.2f} mg·h/L")
    print(f"t_half:   {t_half:.2f} h")

    # 6. Sample concentration-time data
    print("\nConcentration-Time Profile (selected points):")
    print("-" * 30)
    print("Time (h)    Conc (mg/L)")
    print("-" * 30)
    for t in [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]:
        idx = np.where(np.isclose(times, t))[0]
        if len(idx) > 0:
            print(f"  {t:4.0f}        {conc[idx[0]]:.4f}")

    print("\nQuickstart complete!")

    return result


if __name__ == "__main__":
    main()
