#!/usr/bin/env python3
"""
2x2 Crossover Design - Python Example

Run: python python.py
"""

from neopkpd import trial


def main():
    print("2x2 Crossover Study")
    print("=" * 50)

    # Create 2x2 crossover design
    design = trial.crossover_2x2(washout_duration=7.0)

    print(f"Design: {design.n_periods} periods, {design.n_sequences} sequences")
    print(f"Washout: {design.washout_duration} days")
    print(f"Sequences: {design.sequence_assignments}")

    # Create dosing regimens (single dose)
    reference_dose = trial.dosing_single(100.0)  # Single 100 mg reference
    test_dose = trial.dosing_single(100.0)  # Single 100 mg test

    # Create treatment arms (by sequence)
    arms = [
        trial.CrossoverArm("Sequence_AB", 12, treatments=["Reference", "Test"]),
        trial.CrossoverArm("Sequence_BA", 12, treatments=["Test", "Reference"]),
    ]

    # Create trial specification
    spec = trial.TrialSpec(
        name="Bioavailability Crossover Study",
        design=design,
        arms=arms,
        treatments={
            "Reference": reference_dose,
            "Test": test_dose
        },
        pk_sampling_times=[0, 0.25, 0.5, 1, 2, 4, 6, 8, 12, 24],
        endpoints=["cmax", "auc_0_inf"],
        seed=42
    )

    print(f"\nTrial: {spec.name}")
    print("-" * 50)

    # Run simulation
    result = trial.simulate_trial(spec)

    # Display results
    print("\nSequence Results:")
    print("-" * 50)
    for seq_name, seq_result in result.sequences.items():
        pct = seq_result.n_completed / seq_result.n_enrolled * 100
        print(f"  {seq_name}:")
        print(f"    Enrolled:       {seq_result.n_enrolled}")
        print(f"    Completed both: {seq_result.n_completed} ({pct:.1f}%)")

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Overall completion: {result.overall_completion_rate:.1%}")
    print(f"Subjects completing both periods: {result.total_completed}")

    # Period-wise summary
    print("\nPeriod Results:")
    for period in [1, 2]:
        period_result = result.get_period_result(period)
        print(f"  Period {period}: {period_result.n_completed} completed")

    return result


if __name__ == "__main__":
    main()
