#!/usr/bin/env python3
"""
3+3 Dose Escalation - Python Example

Run: python python.py
"""

from neopkpd import trial


def main():
    print("3+3 Dose Escalation Study")
    print("=" * 50)

    # Define dose levels
    dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]

    # Create 3+3 design
    design = trial.dose_escalation_3plus3(
        dose_levels,
        starting_dose=10.0,
        max_dlt_rate=0.33,
        cohort_size=3,
        max_subjects=30
    )

    print(f"Dose levels: {dose_levels} mg")
    print(f"Starting dose: {design.starting_dose} mg")
    print(f"Cohort size: {design.cohort_size}")
    print(f"Max DLT rate: {design.max_dlt_rate:.0%}")
    print(f"Max subjects: {design.max_subjects}")

    # Define DLT probabilities for simulation (true, unknown in real study)
    dlt_probs = {
        10.0: 0.05,
        25.0: 0.10,
        50.0: 0.18,
        100.0: 0.35,
        200.0: 0.55
    }

    # Create trial specification
    spec = trial.EscalationTrialSpec(
        name="Phase I FIH Dose Escalation",
        design=design,
        dlt_probabilities=dlt_probs,  # For simulation only
        regimen_per_dose=lambda dose: trial.dosing_qd(dose, 28),
        seed=42
    )

    print(f"\nTrial: {spec.name}")
    print("-" * 50)

    # Run simulation
    result = trial.simulate_escalation(spec)

    # Display cohort-by-cohort results
    print("\nEscalation History:")
    print("-" * 50)
    for i, cohort in enumerate(result.cohorts, 1):
        decision = cohort.decision
        print(f"Cohort {i} ({cohort.dose:.0f} mg): "
              f"{cohort.n_dlt}/{cohort.n} DLTs -> {decision}")

    print("\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)

    if result.mtd is not None:
        print(f"Recommended MTD: {result.mtd:.0f} mg")
    else:
        print("MTD not determined (study stopped or all doses explored)")

    print(f"Total subjects enrolled: {result.total_subjects}")
    print(f"Highest dose tested: {result.highest_dose_tested:.0f} mg")

    # Summary statistics
    print("\nDose-Level Summary:")
    print("-" * 50)
    print("Dose (mg)  N    DLTs   DLT Rate")
    print("-" * 50)
    for dose_summary in result.dose_summaries:
        rate = dose_summary.n_dlt / dose_summary.n * 100 if dose_summary.n > 0 else 0
        print(f"  {dose_summary.dose:5.0f}    {dose_summary.n:2d}    "
              f"{dose_summary.n_dlt:2d}     {rate:5.1f}%")

    return result


if __name__ == "__main__":
    main()
