#!/usr/bin/env python3
"""
Bioequivalence Study - Python Example

Run: python python.py
"""

from neopkpd import trial
import numpy as np


def main():
    print("Bioequivalence Study")
    print("=" * 50)

    # Create BE design
    design = trial.bioequivalence_design(
        n_periods=2,
        n_sequences=2,
        washout_duration=7.0,
        bioequivalence_limits=(0.80, 1.25),
        parameters=["cmax", "auc_0_inf"],
        regulatory_guidance="fda"
    )

    print(f"Design: {design.n_periods}x{design.n_sequences} crossover")
    print(f"Washout: {design.washout_duration} days")
    print(f"BE limits: {design.bioequivalence_limits}")
    print(f"Parameters: {design.parameters}")

    # Create trial specification
    spec = trial.BETrialSpec(
        name="Generic Bioequivalence Study",
        design=design,
        n_per_sequence=12,  # 24 total
        treatments={
            "Reference": trial.dosing_single(100.0),
            "Test": trial.dosing_single(100.0)
        },
        pk_sampling_times=[0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24],
        seed=42
    )

    print(f"\nTrial: {spec.name}")
    print(f"Total N: {spec.n_per_sequence * 2}")
    print("-" * 50)

    # Simulate BE study with known variability
    # In practice, this would come from the PK model simulation
    result = trial.simulate_be_trial(spec, intra_subject_cv=0.25)

    # Display individual NCA results
    print("\nIndividual Results (first 5 subjects):")
    print("-" * 50)
    print("Subject  Cmax_R   Cmax_T   AUC_R    AUC_T")
    print("-" * 50)
    for i, subj in enumerate(result.subjects[:5]):
        print(f"  {subj.id:3d}    {subj.cmax_ref:.2f}   {subj.cmax_test:.2f}   "
              f"{subj.auc_ref:.1f}   {subj.auc_test:.1f}")

    # Statistical analysis
    print("\n" + "=" * 50)
    print("BIOEQUIVALENCE ASSESSMENT")
    print("=" * 50)

    # Cmax analysis
    cmax_result = trial.assess_bioequivalence(
        test=[s.cmax_test for s in result.subjects],
        reference=[s.cmax_ref for s in result.subjects]
    )

    # AUC analysis
    auc_result = trial.assess_bioequivalence(
        test=[s.auc_test for s in result.subjects],
        reference=[s.auc_ref for s in result.subjects]
    )

    print("\nParameter   GMR      90% CI           Within BE Limits?")
    print("-" * 60)

    cmax_be = "Yes" if cmax_result["bioequivalent"] else "No"
    print(f"Cmax        {cmax_result['point_estimate']:.4f}   "
          f"[{cmax_result['ci_90_lower']:.4f}, {cmax_result['ci_90_upper']:.4f}]   {cmax_be}")

    auc_be = "Yes" if auc_result["bioequivalent"] else "No"
    print(f"AUC0-inf    {auc_result['point_estimate']:.4f}   "
          f"[{auc_result['ci_90_lower']:.4f}, {auc_result['ci_90_upper']:.4f}]   {auc_be}")

    # Overall conclusion
    print("\n" + "=" * 50)
    overall_be = cmax_result["bioequivalent"] and auc_result["bioequivalent"]
    if overall_be:
        print("CONCLUSION: Bioequivalence DEMONSTRATED")
    else:
        print("CONCLUSION: Bioequivalence NOT demonstrated")
    print("=" * 50)

    return {
        "cmax": cmax_result,
        "auc": auc_result,
        "bioequivalent": overall_be
    }


if __name__ == "__main__":
    main()
