#!/usr/bin/env python3
"""
FOCE-I One-Compartment Estimation - Python Example

Run: python python.py
"""

import os
from neopkpd import estimate, create_model_spec, create_observed_data
from neopkpd.estimation import EstimationConfig, FOCEI
import pandas as pd
import numpy as np


def main():
    print("FOCE-I One-Compartment Estimation")
    print("=" * 50)

    # Load observed data
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(script_dir, "data.csv")
    df = pd.read_csv(data_path)

    print(f"\nData Summary:")
    print(f"  Subjects: {df['ID'].nunique()}")
    print(f"  Observations: {len(df)}")
    print(f"  Time range: {df['TIME'].min()} - {df['TIME'].max()} h")

    # Create observed data object
    observed = create_observed_data(
        df,
        id_col="ID",
        time_col="TIME",
        dv_col="DV",
        amt_col="AMT",
        evid_col="EVID"
    )

    # Model specification
    model_spec = create_model_spec(
        "OneCompOralFirstOrder",
        name="foce_example"
    )

    # Estimation configuration
    config = EstimationConfig(
        method=FOCEI(),
        theta_names=["CL", "V", "Ka"],
        theta_init=[4.0, 40.0, 1.0],
        theta_lower=[0.1, 1.0, 0.1],
        theta_upper=[50.0, 500.0, 10.0],
        omega_init={"CL": 0.1, "V": 0.1, "Ka": 0.2},
        omega_structure="diagonal",
        error_model="proportional",
        sigma_init=[0.05],
        maxiter=500,
        abstol=1e-6,
        reltol=1e-4
    )

    # Run estimation
    print("\nRunning FOCE-I estimation...")
    result = estimate(observed, model_spec, config)

    # Display results
    print("\n" + "=" * 60)
    print("ESTIMATION RESULTS")
    print("=" * 60)

    print("\nFixed Effects (θ):")
    print("-" * 50)
    print("Parameter   Estimate    SE        RSE%")
    print("-" * 50)
    for name, est, se in zip(config.theta_names, result["theta"], result["theta_se"]):
        rse = se / est * 100
        print(f"{name:12s}{est:8.3f}    {se:6.3f}    {rse:.1f}%")

    print("\nRandom Effects (ω):")
    print("-" * 50)
    for name, omega in result["omega"].items():
        print(f"ω_{name} = {np.sqrt(omega):.3f} ({np.sqrt(omega)*100:.1f}% CV)")

    print("\nResidual Error (σ):")
    print("-" * 50)
    print(f"σ_prop = {result['sigma'][0]:.4f} ({result['sigma'][0]*100:.1f}% CV)")

    print("\nObjective Function:")
    print("-" * 50)
    print(f"OFV = {result['ofv']:.2f}")
    print(f"AIC = {result['aic']:.2f}")
    print(f"BIC = {result['bic']:.2f}")

    print("\nConvergence:")
    print("-" * 50)
    print(f"Status: {'Converged' if result['converged'] else 'Not converged'}")
    print(f"Iterations: {result['iterations']}")

    # Compare to true values
    print("\n" + "=" * 60)
    print("COMPARISON TO TRUE VALUES")
    print("=" * 60)
    print("Parameter   Estimate   True    Bias%")
    print("-" * 50)
    true_values = [5.0, 50.0, 1.5]
    for name, est, true_val in zip(config.theta_names, result["theta"], true_values):
        bias = (est - true_val) / true_val * 100
        print(f"{name:12s}{est:8.2f}    {true_val}   {bias:.1f}%")

    return result


if __name__ == "__main__":
    main()
