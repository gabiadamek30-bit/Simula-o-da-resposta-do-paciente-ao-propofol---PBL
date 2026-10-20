#!/usr/bin/env python3
"""
Estimation Diagnostic Plots - Python Example

Run: python plot.py
"""

from neopkpd.viz import (
    plot_dv_vs_pred,
    plot_dv_vs_ipred,
    plot_cwres_vs_time,
    plot_cwres_vs_pred,
    plot_eta_distribution,
    plot_gof_panel
)
import matplotlib.pyplot as plt
import numpy as np


def main():
    print("Estimation Diagnostic Plots")
    print("=" * 50)

    # Simulated estimation results (would come from actual estimation)
    np.random.seed(42)
    n_obs = 100

    # Simulated diagnostic data
    dv = np.random.lognormal(2.0, 0.3, n_obs)
    pred = dv * np.random.normal(1.0, 0.1, n_obs)
    ipred = dv * np.random.normal(1.0, 0.05, n_obs)
    time = np.tile(np.array([0, 1, 2, 4, 8, 12, 24, 48, 72, 96]), 10)
    cwres = np.random.normal(0, 1, n_obs)
    etas = {
        "CL": np.random.normal(0, 0.3, 10),
        "V": np.random.normal(0, 0.2, 10),
        "Ka": np.random.normal(0, 0.4, 10)
    }

    # DV vs PRED
    fig, ax = plot_dv_vs_pred(dv, pred, title="DV vs PRED")
    fig.savefig("dv_vs_pred.png", dpi=150, bbox_inches="tight")
    print("Saved: dv_vs_pred.png")

    # DV vs IPRED
    fig, ax = plot_dv_vs_ipred(dv, ipred, title="DV vs IPRED")
    fig.savefig("dv_vs_ipred.png", dpi=150, bbox_inches="tight")
    print("Saved: dv_vs_ipred.png")

    # CWRES vs Time
    fig, ax = plot_cwres_vs_time(cwres, time, title="CWRES vs Time")
    fig.savefig("cwres_vs_time.png", dpi=150, bbox_inches="tight")
    print("Saved: cwres_vs_time.png")

    # CWRES vs PRED
    fig, ax = plot_cwres_vs_pred(cwres, pred, title="CWRES vs PRED")
    fig.savefig("cwres_vs_pred.png", dpi=150, bbox_inches="tight")
    print("Saved: cwres_vs_pred.png")

    # ETA distributions
    fig, axes = plot_eta_distribution(etas, title="ETA Distributions")
    fig.savefig("eta_dist.png", dpi=150, bbox_inches="tight")
    print("Saved: eta_dist.png")

    # Combined GOF panel
    fig = plot_gof_panel(
        dv=dv,
        pred=pred,
        ipred=ipred,
        cwres=cwres,
        time=time,
        title="Goodness-of-Fit"
    )
    fig.savefig("gof_panel.png", dpi=150, bbox_inches="tight")
    print("Saved: gof_panel.png")

    plt.close("all")
    print("\nDone!")


if __name__ == "__main__":
    main()
