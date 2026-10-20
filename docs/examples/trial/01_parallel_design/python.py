#!/usr/bin/env python3
"""
Parallel Group Design - Python Example

Run: python python.py
"""

from neopkpd import trial


def main():
    print("Phase 2 Parallel Study")
    print("=" * 50)

    # Create parallel design
    design = trial.parallel_design(2)  # 2-arm, 1:1 randomization

    print(f"Design: {design.n_arms}-arm parallel")
    print(f"Randomization: {design.randomization_ratio}")

    # Create dosing regimens
    placebo_regimen = trial.dosing_qd(0.0, 28)  # Placebo, 28 days
    active_regimen = trial.dosing_qd(100.0, 28)  # 100 mg QD, 28 days

    # Create treatment arms
    arms = [
        trial.TreatmentArm("Placebo", placebo_regimen, 50, placebo=True),
        trial.TreatmentArm("Active", active_regimen, 50),
    ]

    # Create trial specification
    spec = trial.TrialSpec(
        name="Phase 2 Efficacy Study",
        design=design,
        arms=arms,
        duration_days=28,
        enrollment_rate=5.0,  # 5 subjects per day
        dropout=trial.DropoutSpec(random_rate_per_day=0.005),  # ~14% dropout over 28 days
        compliance=trial.ComplianceSpec(mean_compliance=0.90, pattern="decay"),
        pk_sampling_times=[0, 1, 2, 4, 8, 12, 24],
        endpoints=["pk_exposure"],
        seed=42
    )

    print(f"\nTrial: {spec.name}")
    print(f"Duration: {spec.duration_days} days")
    print("-" * 50)

    # Run simulation
    result = trial.simulate_trial(spec)

    # Display results
    print("\nArm Results:")
    print("-" * 50)
    for arm_name, arm_result in result.arms.items():
        pct = arm_result.n_completed / arm_result.n_enrolled * 100
        print(f"  {arm_name}:")
        print(f"    Enrolled:   {arm_result.n_enrolled}")
        print(f"    Completed:  {arm_result.n_completed} ({pct:.1f}%)")
        print(f"    Dropouts:   {arm_result.n_dropout}")
        print(f"    Compliance: {arm_result.mean_compliance:.1%}")

    print("\n" + "=" * 50)
    print("SUMMARY")
    print("=" * 50)
    print(f"Overall completion: {result.overall_completion_rate:.1%}")
    print(f"Overall compliance: {result.overall_compliance:.1%}")
    print(f"Total enrolled:     {result.total_enrolled}")
    print(f"Total completed:    {result.total_completed}")

    return result


if __name__ == "__main__":
    main()
