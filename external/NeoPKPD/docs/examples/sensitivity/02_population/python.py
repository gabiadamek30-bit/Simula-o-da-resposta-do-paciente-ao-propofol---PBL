#!/usr/bin/env python3
"""
Population Sensitivity - Python Example

Run: python python.py
"""

from neopkpd import create_model_spec, create_population_spec, compute_population_sensitivity


def main():
    print("Population Sensitivity Analysis")
    print("=" * 50)

    # Base model
    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"Ka": 0.16, "CL": 0.09, "V": 0.04}, "seed": 42},
        n=100
    )

    # Compute population sensitivity
    result = compute_population_sensitivity(
        pop_spec,
        parameters=["Ka", "CL", "V", "omega_Ka", "omega_CL", "omega_V"],
        perturbation=0.10,
        metrics=["auc_median", "auc_cv", "cmax_median", "cmax_cv"],
        n_simulations=100,
        t_end=24.0,
        seed=42
    )

    print("\nPopulation Parameters:")
    print("-" * 40)
    print("Fixed Effects (Typical Values):")
    print(f"  Ka = {base_model.params['Ka']} 1/h")
    print(f"  CL = {base_model.params['CL']} L/h")
    print(f"  V  = {base_model.params['V']} L")
    print("\nRandom Effects (Omegas):")
    print(f"  omega_Ka = 0.16 (CV ~40%)")
    print(f"  omega_CL = 0.09 (CV ~30%)")
    print(f"  omega_V  = 0.04 (CV ~20%)")

    print("\nNominal Population Metrics:")
    print("-" * 40)
    for metric, value in result.nominal_metrics.items():
        print(f"  {metric:12s}: {value:.3f}")

    print("\nSensitivity to Typical Values:")
    print("-" * 50)
    print("Parameter    AUC_med   AUC_CV   Cmax_med  Cmax_CV")
    print("-" * 50)
    for param in ["Ka", "CL", "V"]:
        s = result.sensitivity[param]
        print(f"  {param:4s}       {s['auc_median']:+.3f}    {s['auc_cv']:+.3f}    "
              f"{s['cmax_median']:+.3f}    {s['cmax_cv']:+.3f}")

    print("\nSensitivity to IIV (Omegas):")
    print("-" * 50)
    for param in ["omega_Ka", "omega_CL", "omega_V"]:
        s = result.sensitivity[param]
        print(f"  {param:8s}   {s['auc_median']:+.3f}    {s['auc_cv']:+.3f}    "
              f"{s['cmax_median']:+.3f}    {s['cmax_cv']:+.3f}")

    return result


if __name__ == "__main__":
    main()
