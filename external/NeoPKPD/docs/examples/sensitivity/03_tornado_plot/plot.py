#!/usr/bin/env python3
"""
Tornado Plot - Python Example

Run: python plot.py
"""

from neopkpd import create_model_spec, compute_sensitivity
from neopkpd.viz import plot_tornado
import matplotlib.pyplot as plt


def main():
    print("Tornado Plot for Sensitivity Analysis")
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
        perturbation=0.20,  # ±20%
        metrics=["auc", "cmax", "t_half"],
        t_end=24.0
    )

    # Tornado plot for AUC
    fig, ax = plot_tornado(
        result,
        metric="auc",
        title="Sensitivity of AUC to ±20% Parameter Change",
        xlabel="% Change in AUC",
        color_positive="steelblue",
        color_negative="coral"
    )
    fig.savefig("tornado_auc.png", dpi=150, bbox_inches="tight")
    print("Saved: tornado_auc.png")

    # Tornado plot for Cmax
    fig, ax = plot_tornado(
        result,
        metric="cmax",
        title="Sensitivity of Cmax to ±20% Parameter Change",
        xlabel="% Change in Cmax"
    )
    fig.savefig("tornado_cmax.png", dpi=150, bbox_inches="tight")
    print("Saved: tornado_cmax.png")

    # Combined tornado for multiple metrics
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for i, metric in enumerate(["auc", "cmax", "t_half"]):
        plot_tornado(result, metric=metric, ax=axes[i], title=metric.upper())

    fig.tight_layout()
    fig.savefig("tornado_combined.png", dpi=150, bbox_inches="tight")
    print("Saved: tornado_combined.png")

    plt.close("all")
    print("\nDone!")


if __name__ == "__main__":
    main()
