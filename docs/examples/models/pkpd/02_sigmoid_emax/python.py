#!/usr/bin/env python3
"""
Sigmoid Emax (Hill) PD Model - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Sigmoid Emax (Hill) PD Model")
    print("=" * 50)

    # PK parameters
    CL = 5.0   # L/h
    V = 50.0   # L
    Dose = 100.0  # mg

    # PD parameters
    E0 = 0.0      # Baseline effect
    Emax = 100.0  # Maximum effect
    EC50 = 1.0    # mg/L
    gamma = 2.0   # Hill coefficient

    # Create PKPD model specification
    model = create_model_spec(
        "SigmoidEmax",
        name="sigmoid_emax_example",
        params={
            "CL": CL, "V": V,
            "E0": E0, "Emax": Emax, "EC50": EC50, "gamma": gamma
        },
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation
    print("\nRunning PKPD simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=48.0,
        saveat=0.25
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])
    effect = np.array(result["observations"]["effect"])

    # Compute metrics
    cmax = np.max(conc)
    emax_observed = np.max(effect)

    print("\nResults:")
    print("-" * 50)
    print(f"PK Parameters: CL={CL} L/h, V={V} L")
    print("PD Parameters:")
    print(f"  E0 = {E0}, Emax = {Emax}, EC50 = {EC50} mg/L")
    print(f"  gamma (Hill coefficient) = {gamma}")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Max Effect = {emax_observed:.2f}")

    # Show steepness effect
    print(f"\nEffect of Hill Coefficient (gamma = {gamma}):")
    print(f"  At C = 0.5 × EC50: Effect = {E0 + Emax * (0.5**gamma) / (1 + 0.5**gamma):.1f}%")
    print(f"  At C = EC50:       Effect = {E0 + Emax * 0.5:.1f}%")
    print(f"  At C = 2 × EC50:   Effect = {E0 + Emax * (2**gamma) / (1 + 2**gamma):.1f}%")

    # Concentration-Effect profile
    print("\nConcentration-Effect Profile:")
    print("-" * 45)
    print("Time (h)  Conc (mg/L)  Effect")
    print("-" * 45)
    for t in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]:
        idx = np.argmin(np.abs(times - t))
        print(f"  {t:4.0f}     {conc[idx]:7.3f}    {effect[idx]:.2f}")

    return result


if __name__ == "__main__":
    main()
