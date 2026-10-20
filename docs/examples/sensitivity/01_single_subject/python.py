#!/usr/bin/env python3
"""
Single Subject Sensitivity - Python Example

Run: python python.py
"""

from neopkpd import simulate, create_model_spec, compute_sensitivity


def main():
    print("Single Subject Sensitivity Analysis")
    print("=" * 50)

    # Base model
    model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    # Compute sensitivity
    result = compute_sensitivity(
        model,
        parameters=["Ka", "CL", "V"],
        perturbation=0.10,  # ±10%
        metrics=["auc", "cmax", "tmax", "t_half"],
        t_end=24.0
    )

    print("\nNominal Values:")
    print("-" * 40)
    print(f"  Ka = {model.params['Ka']} 1/h")
    print(f"  CL = {model.params['CL']} L/h")
    print(f"  V  = {model.params['V']} L")

    print("\nNominal Metrics:")
    print("-" * 40)
    for metric, value in result.nominal_metrics.items():
        print(f"  {metric:8s}: {value:.3f}")

    print("\nSensitivity Coefficients:")
    print("-" * 50)
    print("Parameter    AUC      Cmax     Tmax     t_half")
    print("-" * 50)
    for param in ["Ka", "CL", "V"]:
        s = result.sensitivity[param]
        print(f"  {param:4s}      {s['auc']:+.3f}   {s['cmax']:+.3f}   "
              f"{s['tmax']:+.3f}   {s['t_half']:+.3f}")

    print("\nInterpretation:")
    print("-" * 50)
    print("Coefficient = (% change in metric) / (% change in parameter)")
    print("  > 0: Positive relationship")
    print("  < 0: Negative relationship")
    print("  ≈ 1: Proportional relationship")

    # Find most influential parameter for each metric
    print("\nMost Influential Parameters:")
    print("-" * 50)
    for metric in ["auc", "cmax", "t_half"]:
        max_param = max(
            ["Ka", "CL", "V"],
            key=lambda p: abs(result.sensitivity[p][metric])
        )
        val = result.sensitivity[max_param][metric]
        print(f"  {metric:8s}: {max_param} (sensitivity = {val:+.3f})")

    return result


if __name__ == "__main__":
    main()
