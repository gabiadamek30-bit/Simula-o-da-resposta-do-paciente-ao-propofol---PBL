#!/usr/bin/env python3
"""
Standard VPC - Python Example

Run: python python.py
"""

from neopkpd import compute_vpc, create_model_spec, create_population_spec, create_observed_data
import numpy as np


def main():
    print("Standard Visual Predictive Check (VPC)")
    print("=" * 50)

    # Create observed data (theophylline-like)
    n_subjects = 12
    times_per_subject = [0, 0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24]

    observed = create_observed_data(
        ids=list(np.repeat(range(1, n_subjects+1), len(times_per_subject))),
        times=times_per_subject * n_subjects,
        dv=[
            10.5, 9.2, 8.5, 7.3, 5.8, 4.2, 3.1, 2.2, 1.5, 0.9, 0.2,
            8.2, 7.5, 7.0, 6.1, 4.9, 3.6, 2.7, 1.9, 1.3, 0.8, 0.15,
            11.8, 10.5, 9.7, 8.3, 6.6, 4.8, 3.5, 2.5, 1.7, 1.0, 0.22,
            9.5, 8.6, 8.0, 6.9, 5.5, 4.0, 3.0, 2.1, 1.4, 0.85, 0.18,
            12.2, 10.8, 10.0, 8.6, 6.9, 5.0, 3.7, 2.6, 1.8, 1.1, 0.25,
            7.8, 7.1, 6.6, 5.7, 4.6, 3.3, 2.5, 1.8, 1.2, 0.72, 0.13,
            10.1, 9.1, 8.4, 7.2, 5.7, 4.1, 3.1, 2.2, 1.5, 0.88, 0.19,
            9.8, 8.9, 8.2, 7.1, 5.6, 4.1, 3.0, 2.1, 1.45, 0.87, 0.18,
            11.2, 10.0, 9.2, 7.9, 6.3, 4.6, 3.4, 2.4, 1.6, 0.98, 0.21,
            8.8, 8.0, 7.4, 6.4, 5.1, 3.7, 2.8, 1.95, 1.32, 0.79, 0.16,
            10.8, 9.7, 9.0, 7.7, 6.1, 4.4, 3.3, 2.3, 1.55, 0.94, 0.20,
            9.2, 8.3, 7.7, 6.6, 5.3, 3.8, 2.85, 2.0, 1.38, 0.82, 0.17
        ]
    )

    # Population model
    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 0.5, "V": 30.0},
        doses=[{"time": 0.0, "amount": 320.0}]
    )

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"Ka": 0.4, "CL": 0.25, "V": 0.15}, "seed": 12345},
        n=12
    )

    # Generate VPC
    print("\nGenerating VPC (500 simulations)...")
    vpc_result = compute_vpc(
        observed,
        pop_spec,
        n_simulations=500,
        quantiles=[0.05, 0.5, 0.95],
        seed=12345
    )

    print("\n" + "=" * 60)
    print("VPC RESULTS")
    print("=" * 60)

    print("\nObserved Data Quantiles:")
    print("-" * 50)
    print("Time (h)  5th%   Median  95th%")
    print("-" * 50)
    for i, t in enumerate(vpc_result["time_bins"]):
        obs_5 = vpc_result["observed_quantiles"]["0.05"][i]
        obs_50 = vpc_result["observed_quantiles"]["0.5"][i]
        obs_95 = vpc_result["observed_quantiles"]["0.95"][i]
        print(f"  {t:4.1f}   {obs_5:.2f}   {obs_50:.2f}    {obs_95:.2f}")

    pct_within = sum(vpc_result["obs_within_sim_ci"]) / len(vpc_result["obs_within_sim_ci"]) * 100
    print(f"\n% Observations within simulated CI: {pct_within:.1f}%")

    if pct_within > 80:
        print("\nConclusion: VPC indicates adequate model fit")
    else:
        print("\nConclusion: VPC suggests model misspecification")

    return vpc_result


if __name__ == "__main__":
    main()
