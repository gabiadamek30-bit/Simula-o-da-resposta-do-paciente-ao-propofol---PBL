#!/usr/bin/env python3
"""
Time-Varying Covariates (Renal Function) - Python Example

Run: python python.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
import numpy as np


def main():
    print("Time-Varying Covariates (Renal Function)")
    print("=" * 50)

    np.random.seed(12345)

    n_subjects = 50

    def generate_crcl_trajectory(baseline_crcl, decline_rate, times):
        return [max(baseline_crcl - decline_rate * t, 20.0) for t in times]

    covariate_times = list(np.arange(0.0, 169.0, 24.0))

    covariates_data = []
    for i in range(n_subjects):
        baseline_crcl = 100.0 + 20.0 * np.random.randn()
        decline_rate = max(0, 7.0 + 3.0 * np.random.randn())
        crcl_values = generate_crcl_trajectory(baseline_crcl, decline_rate, covariate_times)
        covariates_data.append({
            "CRCL": {str(t): v for t, v in zip(covariate_times, crcl_values)}
        })

    base_model = create_model_spec(
        "OneCompIVBolus",
        name="time_varying_cov_example",
        params={"CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    covariate_model = {
        "effects": [
            {"covariate": "CRCL", "parameter": "CL", "type": "power",
             "exponent": 0.5, "reference": 100.0, "time_varying": True}
        ]
    }

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"CL": 0.2, "V": 0.15}, "seed": 12345},
        covariate_model=covariate_model,
        covariates=covariates_data,
        n=n_subjects
    )

    print(f"\nSimulating {n_subjects} subjects over 7 days...")
    print("CrCL declines over time, affecting drug clearance")

    result = simulate_population(pop_spec, t_start=0.0, t_end=168.0, saveat=2.0)

    print("\nExample Subject CrCL Trajectory:")
    print("-" * 40)
    print("Day   CrCL (mL/min)  CL adjustment")
    print("-" * 40)
    for t in [0.0, 48.0, 96.0, 168.0]:
        crcl = covariates_data[0]["CRCL"].get(str(t), 50.0)
        cl_adj = (crcl / 100.0)**0.5
        print(f"  {int(t/24)}     {crcl:.1f}         {cl_adj:.3f}")

    print("\nNote: Drug accumulation may increase as CrCL decreases")

    return result


if __name__ == "__main__":
    main()
