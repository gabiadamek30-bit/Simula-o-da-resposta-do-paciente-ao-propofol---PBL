#!/usr/bin/env python3
"""
Biophase Equilibration (Effect Compartment) Model - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Biophase Equilibration (Effect Compartment) Model")
    print("=" * 50)

    # PK parameters
    CL = 5.0   # L/h
    V = 50.0   # L
    Dose = 100.0  # mg

    # PD parameters
    Ke0 = 0.5     # Effect site equilibration rate (1/h)
    E0 = 0.0      # Baseline effect
    Emax = 100.0  # Maximum effect
    EC50 = 1.0    # mg/L (at effect site)

    # Calculate equilibration half-life
    t_half_ke0 = np.log(2) / Ke0

    # Create PKPD model specification
    model = create_model_spec(
        "BiophaseEquilibration",
        name="biophase_equilibration_example",
        params={
            "CL": CL, "V": V,
            "Ke0": Ke0, "E0": E0, "Emax": Emax, "EC50": EC50
        },
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation
    print("\nRunning PKPD simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=48.0,
        saveat=0.1
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])  # Plasma concentration
    ce = np.array(result["observations"]["Ce"])       # Effect site concentration
    effect = np.array(result["observations"]["effect"])

    # Find peaks
    cmax = np.max(conc)
    tmax_pk = times[np.argmax(conc)]
    emax_observed = np.max(effect)
    tmax_pd = times[np.argmax(effect)]
    delay = tmax_pd - tmax_pk

    print("\nResults:")
    print("-" * 50)
    print(f"PK Parameters: CL={CL} L/h, V={V} L")
    print("PD Parameters:")
    print(f"  Ke0 = {Ke0} 1/h (t½ke0 = {t_half_ke0:.2f} h)")
    print(f"  E0 = {E0}, Emax = {Emax}, EC50 = {EC50} mg/L")
    print("\nPK:")
    print(f"  Cmax (plasma) = {cmax:.4f} mg/L at t = {tmax_pk:.2f} h")
    print("\nPD:")
    print(f"  Max Effect = {emax_observed:.2f} at t = {tmax_pd:.2f} h")
    print(f"  Effect delay = {delay:.2f} h")

    # Show hysteresis
    print("\nPlasma vs Effect Site Concentration (hysteresis):")
    print("-" * 55)
    print("Time (h)  Cp (mg/L)  Ce (mg/L)  Effect")
    print("-" * 55)
    for t in [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]:
        idx = np.argmin(np.abs(times - t))
        print(f"  {t:4.0f}     {conc[idx]:7.3f}   {ce[idx]:7.3f}    {effect[idx]:.1f}")

    print("\nNote: Counter-clockwise hysteresis - Ce lags behind Cp")
    print("At early times: Cp > Ce (effect building)")
    print("At late times: Ce > Cp (effect persisting)")

    return result


if __name__ == "__main__":
    main()
