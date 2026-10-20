#!/usr/bin/env python3
"""
Basic NCA Analysis - Python Example

Run: python python.py
"""

from neopkpd.nca import compute_nca
import numpy as np


def main():
    print("Basic Non-Compartmental Analysis (NCA)")
    print("=" * 50)

    # Sample concentration-time data (100 mg IV bolus)
    times = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
    concentrations = [2.0, 1.81, 1.64, 1.35, 0.91, 0.61, 0.41, 0.18, 0.03]
    dose = 100.0  # mg

    print("\nConcentration-Time Data:")
    print("-" * 30)
    print("Time (h)   Conc (mg/L)")
    print("-" * 30)
    for t, c in zip(times, concentrations):
        print(f"  {t:4.1f}       {c:.3f}")

    # Compute NCA metrics
    result = compute_nca(
        times=times,
        concentrations=concentrations,
        dose=dose,
        route="iv",
        method="linear_log"
    )

    print("\n" + "=" * 50)
    print("NCA RESULTS")
    print("=" * 50)

    print("\nExposure Metrics:")
    print("-" * 40)
    print(f"Cmax      = {result['cmax']:.3f} mg/L")
    print(f"Tmax      = {result['tmax']:.2f} h")
    print(f"Clast     = {result['clast']:.4f} mg/L")
    print(f"Tlast     = {result['tlast']:.1f} h")

    print("\nAUC:")
    print("-" * 40)
    print(f"AUC_0_t   = {result['auc_0_t']:.2f} mg·h/L")
    print(f"AUC_0_inf = {result['auc_0_inf']:.2f} mg·h/L")
    print(f"AUC_%ext  = {result['auc_pct_extrapolated']:.1f}%")

    print("\nElimination:")
    print("-" * 40)
    print(f"λz        = {result['lambda_z']:.4f} 1/h")
    print(f"t_half    = {result['t_half']:.2f} h")
    print(f"R²        = {result['r_squared']:.4f}")
    print(f"N points  = {result['n_points_lambda_z']}")

    print("\nClearance and Volume:")
    print("-" * 40)
    print(f"CL        = {result['cl']:.2f} L/h")
    print(f"Vz        = {result['vz']:.1f} L")
    print(f"Vss       = {result['vss']:.1f} L")

    print("\nMean Residence Time:")
    print("-" * 40)
    print(f"MRT       = {result['mrt']:.2f} h")
    print(f"AUMC_0_inf= {result['aumc_0_inf']:.1f} mg·h²/L")

    return result


if __name__ == "__main__":
    main()
