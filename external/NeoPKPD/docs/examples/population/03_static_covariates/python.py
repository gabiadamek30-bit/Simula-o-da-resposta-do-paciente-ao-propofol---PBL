#!/usr/bin/env python3
"""
Static Covariates (Weight, Age) - Python Example

Run: python python.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
import numpy as np


def main():
    print("Static Covariates (Weight, Age)")
    print("=" * 50)

    np.random.seed(12345)

    # Generate covariate data
    n_subjects = 100
    weights = np.clip(70.0 + 15.0 * np.random.randn(n_subjects), 40.0, 120.0)
    ages = np.clip(50.0 + 15.0 * np.random.randn(n_subjects), 18.0, 80.0)

    # Base model
    base_model = create_model_spec(
        "OneCompIVBolus",
        name="covariate_example",
        params={"CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    # Covariate model
    covariate_model = {
        "effects": [
            {"covariate": "WT", "parameter": "CL", "type": "power", "exponent": 0.75, "reference": 70.0},
            {"covariate": "WT", "parameter": "V", "type": "power", "exponent": 1.0, "reference": 70.0},
            {"covariate": "AGE", "parameter": "CL", "type": "linear", "slope": -0.01, "reference": 40.0}
        ]
    }

    # Create covariates
    covariates = [{"WT": float(weights[i]), "AGE": float(ages[i])} for i in range(n_subjects)]

    # Population specification
    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"CL": 0.25, "V": 0.15}, "seed": 12345},
        covariate_model=covariate_model,
        covariates=covariates,
        n=n_subjects
    )

    print("\nRunning population simulation with covariates...")
    result = simulate_population(pop_spec, t_start=0.0, t_end=24.0, saveat=0.5)

    print("\nCovariate Distribution:")
    print("-" * 40)
    print(f"Weight: mean={np.mean(weights):.1f} kg, SD={np.std(weights):.1f}")
    print(f"Age:    mean={np.mean(ages):.1f} y, SD={np.std(ages):.1f}")

    individual_params = result["params"]
    CL_values = np.array([p["CL"] for p in individual_params])
    V_values = np.array([p["V"] for p in individual_params])

    print("\nIndividual Parameters (after covariate adjustment):")
    print("-" * 50)
    print(f"CL: mean={np.mean(CL_values):.2f} L/h, range=[{np.min(CL_values):.2f}, {np.max(CL_values):.2f}]")
    print(f"V:  mean={np.mean(V_values):.1f} L, range=[{np.min(V_values):.1f}, {np.max(V_values):.1f}]")

    print("\nExample Covariate Effects:")
    print("-" * 50)
    print(f"50 kg, 30y: CL = 5.0 × (50/70)^0.75 × (1-0.01×(30-40)) = {5.0 * (50/70)**0.75 * (1 - 0.01*(30-40)):.2f} L/h")
    print(f"90 kg, 70y: CL = 5.0 × (90/70)^0.75 × (1-0.01×(70-40)) = {5.0 * (90/70)**0.75 * (1 - 0.01*(70-40)):.2f} L/h")

    return result


if __name__ == "__main__":
    main()
