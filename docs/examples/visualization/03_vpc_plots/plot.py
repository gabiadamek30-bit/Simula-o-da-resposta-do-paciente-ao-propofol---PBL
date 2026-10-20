#!/usr/bin/env python3
"""
VPC Plot - Python Example

Run: python plot.py
"""

from neopkpd import compute_vpc, create_model_spec, create_population_spec, create_observed_data
from neopkpd.viz import plot_vpc
import matplotlib.pyplot as plt
import numpy as np


def main():
    print("VPC Plot")
    print("=" * 50)

    # Create observed data
    n_subjects = 12
    times = [0, 0.25, 0.5, 1, 2, 4, 8, 12, 24]

    observed = create_observed_data(
        ids=list(np.repeat(range(1, n_subjects+1), len(times))),
        times=times * n_subjects,
        dv=np.random.lognormal(mean=2.0, sigma=0.3, size=n_subjects*len(times)).tolist()
    )

    # Population model
    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"Ka": 0.16, "CL": 0.09, "V": 0.04}, "seed": 42},
        n=n_subjects
    )

    # Compute VPC
    vpc_result = compute_vpc(
        observed,
        pop_spec,
        n_simulations=500,
        quantiles=[0.05, 0.5, 0.95],
        seed=42
    )

    # Standard VPC plot
    fig, ax = plot_vpc(
        vpc_result,
        title="Visual Predictive Check",
        obs_marker="o",
        obs_color="black",
        obs_alpha=0.5,
        sim_median_color="blue",
        sim_ci_color="lightblue"
    )
    fig.savefig("vpc.png", dpi=150, bbox_inches="tight")
    print("Saved: vpc.png")

    # Semi-log VPC
    fig, ax = plot_vpc(
        vpc_result,
        title="VPC (Semi-log)",
        log_y=True
    )
    fig.savefig("vpc_semilog.png", dpi=150, bbox_inches="tight")
    print("Saved: vpc_semilog.png")

    plt.close("all")
    print("\nDone!")


if __name__ == "__main__":
    main()
