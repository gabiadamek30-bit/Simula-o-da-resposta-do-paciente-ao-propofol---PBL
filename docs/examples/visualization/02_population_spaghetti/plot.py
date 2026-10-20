#!/usr/bin/env python3
"""
Population Spaghetti Plot - Python Example

Run: python plot.py
"""

from neopkpd import simulate_population, create_model_spec, create_population_spec
from neopkpd.viz import plot_population_spaghetti
import matplotlib.pyplot as plt


def main():
    print("Population Spaghetti Plot")
    print("=" * 50)

    # Create population model
    base_model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    pop_spec = create_population_spec(
        base_model_spec=base_model,
        iiv={"kind": "LogNormalIIV", "omegas": {"Ka": 0.16, "CL": 0.09, "V": 0.04}, "seed": 42},
        n=20
    )

    result = simulate_population(pop_spec, t_end=24.0, saveat=0.5)

    # Spaghetti plot
    fig, ax = plot_population_spaghetti(
        result,
        title="Population PK (n=20)",
        individual_alpha=0.3,
        show_median=True,
        median_color="red",
        median_linewidth=2
    )
    fig.savefig("spaghetti.png", dpi=150, bbox_inches="tight")
    print("Saved: spaghetti.png")

    # With confidence bands
    fig, ax = plot_population_spaghetti(
        result,
        title="Population PK with 90% PI",
        show_ci=True,
        ci_levels=[0.05, 0.95],
        ci_alpha=0.2
    )
    fig.savefig("spaghetti_ci.png", dpi=150, bbox_inches="tight")
    print("Saved: spaghetti_ci.png")

    plt.close("all")
    print("\nDone!")


if __name__ == "__main__":
    main()
