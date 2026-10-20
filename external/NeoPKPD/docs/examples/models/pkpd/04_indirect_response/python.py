#!/usr/bin/env python3
"""
Indirect Response (Turnover) Model - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec
import numpy as np


def main():
    print("Indirect Response (Turnover) Model - Type I")
    print("=" * 50)

    # PK parameters
    CL = 5.0   # L/h
    V = 50.0   # L
    Dose = 100.0  # mg

    # PD parameters (IRM-I: inhibition of production)
    Kin = 10.0    # Production rate (units/h)
    Kout = 0.1    # Loss rate (1/h)
    Imax = 0.9    # Maximum inhibition (90%)
    IC50 = 1.0    # Potency (mg/L)

    # Calculate derived parameters
    R0 = Kin / Kout  # Baseline response
    t_half_R = np.log(2) / Kout  # Response half-life

    # Create PKPD model specification
    model = create_model_spec(
        "IndirectResponseI",  # Type I = inhibition of Kin
        name="indirect_response_example",
        params={
            "CL": CL, "V": V,
            "Kin": Kin, "Kout": Kout, "Imax": Imax, "IC50": IC50
        },
        doses=[{"time": 0.0, "amount": Dose}]
    )

    # Run simulation - 4 days
    print("\nRunning PKPD simulation...")
    result = simulate(
        model,
        t_start=0.0,
        t_end=96.0,
        saveat=0.5
    )

    # Extract results
    times = np.array(result["t"])
    conc = np.array(result["observations"]["conc"])
    response = np.array(result["observations"]["response"])

    # Compute metrics
    cmax = np.max(conc)
    min_response = np.min(response)
    percent_decrease = (R0 - min_response) / R0 * 100
    time_to_min = times[np.argmin(response)]

    print("\nResults:")
    print("-" * 50)
    print(f"PK Parameters: CL={CL} L/h, V={V} L")
    print("PD Parameters (Indirect Response Type I):")
    print(f"  Kin = {Kin} units/h, Kout = {Kout} 1/h")
    print(f"  Baseline R0 = Kin/Kout = {R0:.1f} units")
    print(f"  Response t½ = {t_half_R:.1f} h")
    print(f"  Imax = {Imax} ({int(Imax*100)}% max inhibition)")
    print(f"  IC50 = {IC50} mg/L")
    print("\nSimulated:")
    print(f"  Cmax = {cmax:.4f} mg/L")
    print(f"  Minimum Response = {min_response:.2f} units")
    print(f"  Max % Decrease = {percent_decrease:.1f}%")
    print(f"  Time to nadir = {time_to_min:.1f} h")

    # Show turnover dynamics
    print("\nPK/PD Time Course:")
    print("-" * 55)
    print("Time (h)  Conc (mg/L)  I(C)      Response")
    print("-" * 55)
    for t in [0.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0, 96.0]:
        idx = np.argmin(np.abs(times - t))
        c = conc[idx]
        inhibition = Imax * c / (IC50 + c)
        r = response[idx]
        print(f"  {t:4.0f}      {c:7.3f}  {inhibition:5.2f}     {r:.1f}")

    print("\nNote: Response delayed relative to concentration")
    print(f"- Drug suppresses Kin → Response decreases slowly (t½={t_half_R:.1f}h)")
    print("- After drug washout → Response returns to baseline")

    return result


if __name__ == "__main__":
    main()
