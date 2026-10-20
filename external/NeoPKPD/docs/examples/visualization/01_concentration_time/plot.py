#!/usr/bin/env python3
"""
Concentration-Time Plot - Python Example

Run: python plot.py
"""

from neopkpd import simulate, create_model_spec
from neopkpd.viz import plot_concentration_time
import matplotlib.pyplot as plt


def main():
    print("Concentration-Time Plot")
    print("=" * 50)

    # Create model and simulate
    model = create_model_spec(
        "OneCompOralFirstOrder",
        params={"Ka": 1.5, "CL": 5.0, "V": 50.0},
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    result = simulate(model, t_end=24.0, saveat=0.1)

    # Basic plot
    fig, ax = plot_concentration_time(
        result,
        title="Oral PK Profile (100 mg)",
        xlabel="Time (hours)",
        ylabel="Concentration (ng/mL)"
    )
    fig.savefig("pk_linear.png", dpi=150, bbox_inches="tight")
    print("Saved: pk_linear.png")

    # Semi-log plot
    fig, ax = plot_concentration_time(
        result,
        title="Oral PK Profile (Semi-log)",
        xlabel="Time (hours)",
        ylabel="Concentration (ng/mL)",
        log_y=True
    )
    fig.savefig("pk_semilog.png", dpi=150, bbox_inches="tight")
    print("Saved: pk_semilog.png")

    # With dose annotation
    fig, ax = plot_concentration_time(
        result,
        title="PK Profile with Dose",
        show_doses=True,
        dose_color="red"
    )
    fig.savefig("pk_with_dose.png", dpi=150, bbox_inches="tight")
    print("Saved: pk_with_dose.png")

    plt.close("all")
    print("\nDone!")


if __name__ == "__main__":
    main()
