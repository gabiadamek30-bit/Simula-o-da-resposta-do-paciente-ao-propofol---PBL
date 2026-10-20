#!/usr/bin/env python3
"""
Transit Compartment Absorption Model - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Transit Compartment Absorption Model")
    print("=" * 50)

    # Model parameters
    Ktr = 2.0   # Transit rate constant (1/h)
    n = 3       # Number of transit compartments
    CL = 5.0    # Clearance (L/h)
    V = 50.0    # Volume of distribution (L)
    Dose = 100.0  # mg

    # Calculate mean transit time
    MTT = (n + 1) / Ktr

    # Create model specification
    model = create_model_spec(
        "TransitAbsorption",
        name="transit_absorption_example",
        params={"Ktr": Ktr, "n": n, "CL": CL, "V": V},
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation
    print("\nRunning simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=48.0,
        saveat=0.1
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])

    # Compute metrics
    cmax = np.max(conc)
    tmax = times[np.argmax(conc)]
    auc = np.trapz(conc, times)

    # Compare with first-order absorption
    Ka_equivalent = 1.0 / MTT

    print("\nResults:")
    print("-" * 50)
    print("Transit Model Parameters:")
    print(f"  Ktr = {Ktr} 1/h")
    print(f"  n   = {n} compartments")
    print(f"  MTT = {MTT:.2f} h")
    print("\nPK Parameters:")
    print(f"  CL = {CL} L/h")
    print(f"  V  = {V} L")
    print("\nSimulated PK:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Tmax = {tmax:.2f} h")
    print(f"  AUC  = {auc:.2f} mg·h/L")
    print(f"\nNote: First-order model with equivalent MTT would have Ka ≈ {Ka_equivalent:.2f} 1/h")

    # Show absorption profile
    print("\nConcentration-Time Profile:")
    print("-" * 40)
    for t in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]:
        idx = np.argmin(np.abs(times - t))
        marker = " <-- Tmax" if np.isclose(t, tmax, atol=0.1) else ""
        print(f"t={t:4.1f}h: {conc[idx]:.4f} mg/L{marker}")

    return result


if __name__ == "__main__":
    main()
