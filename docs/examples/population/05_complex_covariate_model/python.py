#!/usr/bin/env python3
"""
Complex Covariate Model - Python Example

Run: python python.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
import numpy as np


def main():
    print("Complex Covariate Model")
    print("=" * 50)

    np.random.seed(12345)

    n_subjects = 200

    # Generate covariates
    weights = np.clip(70.0 + 15.0 * np.random.randn(n_subjects), 40.0, 120.0)
    ages = np.clip(50.0 + 15.0 * np.random.randn(n_subjects), 18.0, 80.0)
    crcl = np.clip(100.0 + 25.0 * np.random.randn(n_subjects), 30.0, 150.0)
    sex = np.random.randint(0, 2, n_subjects)

    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        name="complex_covariate_example",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    covariate_model = {
        "effects": [
            {"covariate": "WT", "parameter": "CL", "type": "power", "exponent": 0.75, "reference": 70.0},
            {"covariate": "WT", "parameter": "V", "type": "power", "exponent": 1.0, "reference": 70.0},
            {"covariate": "CRCL", "parameter": "CL", "type": "power", "exponent": 0.5, "reference": 100.0},
            {"covariate": "SEX", "parameter": "CL", "type": "categorical", "coefficient": 0.85}
        ]
    }

    covariates = [{"WT": float(weights[i]), "AGE": float(ages[i]),
                   "CRCL": float(crcl[i]), "SEX": int(sex[i])} for i in range(n_subjects)]

    omega_matrix = [[0.09, 0.02, 0.0], [0.02, 0.04, 0.0], [0.0, 0.0, 0.16]]

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={
            "kind": "LogNormalIIVCorrelated",
            "parameters": ["CL", "V", "Ka"],
            "omega_matrix": omega_matrix,
            "seed": 12345
        },
        covariate_model=covariate_model,
        covariates=covariates,
        n=n_subjects
    )

    print(f"\nRunning complex population simulation...")
    print(f"  Subjects: {n_subjects}")
    print("  Covariates: WT, AGE, CrCL, SEX")
    print("  Correlated IIV: CL-V correlation")

    result = simulate_population(pop_spec, t_start=0.0, t_end=24.0, saveat=0.5)

    individual_params = result["params"]
    CL_values = np.array([p["CL"] for p in individual_params])
    V_values = np.array([p["V"] for p in individual_params])

    print("\nPopulation Summary:")
    print("-" * 50)
    print("Parameter   Mean    SD      CV%")
    print("-" * 50)
    print(f"CL (L/h)    {np.mean(CL_values):.2f}   {np.std(CL_values):.2f}   {np.std(CL_values)/np.mean(CL_values)*100:.1f}")
    print(f"V (L)       {np.mean(V_values):.1f}   {np.std(V_values):.1f}   {np.std(V_values)/np.mean(V_values)*100:.1f}")

    CL_male = CL_values[sex == 0]
    CL_female = CL_values[sex == 1]
    print("\nCL by Sex:")
    print(f"  Male:   {np.mean(CL_male):.2f} L/h (n={len(CL_male)})")
    print(f"  Female: {np.mean(CL_female):.2f} L/h (n={len(CL_female)})")
    print(f"  Ratio:  {np.mean(CL_female)/np.mean(CL_male):.2f}")

    return result


if __name__ == "__main__":
    main()
