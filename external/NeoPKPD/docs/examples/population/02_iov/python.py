#!/usr/bin/env python3
"""
Inter-Occasion Variability (IOV) - Python Example

Run: python python.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
import numpy as np


def main():
    print("Inter-Occasion Variability (IOV)")
    print("=" * 50)

    # Base model
    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        name="iov_example",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    # Population specification with IIV and IOV
    n_subjects = 50
    n_occasions = 3
    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={
            "kind": "LogNormalIIV",
            "omegas": {"Ka": 0.4, "CL": 0.3, "V": 0.2},
            "seed": 12345
        },
        iov={
            "kind": "LogNormalIOV",
            "pis": {"Ka": 0.2, "CL": 0.15},
            "occasions": n_occasions,
            "seed": 54321
        },
        n=n_subjects
    )

    print(f"\nRunning population simulation...")
    print(f"  Subjects: {n_subjects}")
    print(f"  Occasions: {n_occasions}")

    result = simulate_population(
        pop_spec,
        t_start=0.0,
        t_end=24.0,
        saveat=0.5
    )

    print("\nVariance Components:")
    print("-" * 50)
    print("Parameter   ω² (IIV)   π² (IOV)   Total Var")
    print("-" * 50)
    print("Ka          0.16       0.04       0.20")
    print("CL          0.09       0.0225     0.1125")
    print("V           0.04       -          0.04")

    print("\nConclusion: IOV adds to total variability within subjects")
    print("Useful for modeling replicate PK studies")

    return result


if __name__ == "__main__":
    main()
