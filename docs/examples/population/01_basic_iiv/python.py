#!/usr/bin/env python3
"""
Basic Inter-Individual Variability (IIV) - Python Example

Run: python python.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
import numpy as np


def main():
    print("Basic Inter-Individual Variability (IIV)")
    print("=" * 50)

    # Base PK model
    base_model = create_model_spec(
        "OneCompIVBolus",
        name="population_iiv_example",
        params={"CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    # Population specification with IIV
    n_subjects = 100
    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={
            "kind": "LogNormalIIV",
            "omegas": {"CL": 0.3, "V": 0.2},
            "seed": 12345
        },
        n=n_subjects
    )

    # Run population simulation
    print(f"\nRunning population simulation (n={n_subjects})...")
    result = simulate_population(
        pop_spec,
        t_start=0.0,
        t_end=24.0,
        saveat=0.5
    )

    # Extract individual parameters
    individual_params = result["params"]
    CL_values = np.array([p["CL"] for p in individual_params])
    V_values = np.array([p["V"] for p in individual_params])

    print("\nIndividual Parameter Distribution:")
    print("-" * 50)
    print("Parameter  Mean    Median   SD      5th%    95th%")
    print("-" * 50)
    print(f"CL (L/h)   {np.mean(CL_values):.2f}   {np.median(CL_values):.2f}    {np.std(CL_values):.2f}   {np.percentile(CL_values, 5):.2f}   {np.percentile(CL_values, 95):.2f}")
    print(f"V (L)      {np.mean(V_values):.1f}   {np.median(V_values):.1f}   {np.std(V_values):.1f}   {np.percentile(V_values, 5):.1f}   {np.percentile(V_values, 95):.1f}")

    # Population summary of concentrations
    summaries = result["summaries"]
    conc_summary = summaries["conc"]
    times = np.array(result["t"])

    print("\nPopulation Concentration Summary:")
    print("-" * 50)
    print("Time (h)  Mean    Median  5th%    95th%")
    print("-" * 50)
    for t in [0.0, 2.0, 4.0, 8.0, 12.0, 24.0]:
        idx = np.argmin(np.abs(times - t))
        mean_val = conc_summary["mean"][idx]
        median_val = conc_summary["median"][idx]
        p5 = conc_summary["quantiles"]["0.05"][idx]
        p95 = conc_summary["quantiles"]["0.95"][idx]
        print(f"  {t:4.0f}     {mean_val:.3f}  {median_val:.3f}  {p5:.3f}  {p95:.3f}")

    print("\nNote: IIV creates a distribution of PK profiles")
    print("Each subject has unique CL and V values")

    return result


if __name__ == "__main__":
    main()
