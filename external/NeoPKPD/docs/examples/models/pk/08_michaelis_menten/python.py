#!/usr/bin/env python3
"""
Michaelis-Menten Elimination Model - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Michaelis-Menten Elimination Model")
    print("=" * 50)

    # Model parameters
    Vmax = 50.0  # Maximum elimination rate (mg/h)
    Km = 10.0    # Michaelis constant (mg/L)
    V = 50.0     # Volume of distribution (L)
    Dose = 500.0 # mg (high dose to show saturation)

    # Calculate derived parameters
    CLint = Vmax / Km  # Intrinsic clearance at low concentrations
    t_half_linear = 0.693 * V * Km / Vmax

    # Create model specification
    model = create_model_spec(
        "MichaelisMentenElimination",
        name="michaelis_menten_example",
        params={"Vmax": Vmax, "Km": Km, "V": V},
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation
    print("\nRunning simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=72.0,
        saveat=0.25
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])

    # Compute metrics
    cmax = np.max(conc)
    auc = np.trapz(conc, times)

    # Analyze saturation
    C0 = Dose / V
    print("\nResults:")
    print("-" * 50)
    print("Model Parameters:")
    print(f"  Vmax = {Vmax} mg/h")
    print(f"  Km   = {Km} mg/L")
    print(f"  V    = {V} L")
    print("\nDerived Parameters:")
    print(f"  CLint (at low C) = {CLint:.2f} L/h")
    print(f"  t½ (at low C)    = {t_half_linear:.2f} h")
    print("\nSaturation Analysis:")
    print(f"  C0 = {C0:.2f} mg/L (Dose/V)")
    print(f"  C0/Km = {C0/Km:.1f} (>>1 = saturated)")
    print("\nSimulated PK:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  AUC  = {auc:.2f} mg·h/L")

    # Show elimination kinetics changing
    print("\nConcentration-Time Profile (nonlinear elimination):")
    print("-" * 55)
    print("Time (h)  Conc (mg/L)  C/Km    Kinetics")
    print("-" * 55)
    for t in [0.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0]:
        idx = np.argmin(np.abs(times - t))
        c = conc[idx]
        c_km = c / Km
        if c_km > 5:
            kinetics = "zero-order (saturated)"
        elif c_km > 0.2:
            kinetics = "mixed"
        else:
            kinetics = "first-order"
        print(f"  {t:4.0f}      {c:7.3f}    {c_km:5.2f}   {kinetics}")

    # Compare to linear model
    print(f"\nNote: A linear model with CL = CLint = {CLint:.1f} L/h")
    print(f"would have t½ = {t_half_linear:.1f} h and AUC = {Dose/CLint:.1f} mg·h/L")
    print(f"The actual AUC is {(auc/(Dose/CLint) * 100 - 100):.0f}% higher due to saturation.")

    return result


if __name__ == "__main__":
    main()
