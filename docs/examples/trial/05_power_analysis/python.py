#!/usr/bin/env python3
"""
Power Analysis - Python Example

Run: python python.py
"""

from neopkpd import trial
import math


def main():
    print("Power Analysis and Sample Size Estimation")
    print("=" * 60)

    # 1. Power calculation for given sample size
    print("\n1. POWER CALCULATION")
    print("-" * 40)

    result = trial.estimate_power_analytical(
        n_per_arm=50,
        effect_size=0.5,  # Cohen's d (medium effect)
        sd=1.0,
        alpha=0.05,
        alternative="two-sided"
    )

    print(f"Given:")
    print(f"  N per arm:    50")
    print(f"  Effect size:  0.5 (medium)")
    print(f"  Alpha:        0.05")
    print(f"\nCalculated power: {result.power:.1%}")

    # 2. Sample size for target power
    print("\n" + "=" * 60)
    print("2. SAMPLE SIZE ESTIMATION")
    print("-" * 40)

    ss_result = trial.estimate_sample_size(
        target_power=0.80,
        effect_size=0.5,
        alpha=0.05
    )

    print(f"Target power: 80%")
    print(f"Effect size:  0.5 (medium)")
    print(f"Alpha:        0.05")
    print(f"\nRequired N per arm: {ss_result.n_per_arm}")
    print(f"Total N:            {ss_result.total_n}")
    print(f"Achieved power:     {ss_result.achieved_power:.1%}")

    # 3. Adjust for dropout
    print("\n" + "=" * 60)
    print("3. DROPOUT ADJUSTMENT")
    print("-" * 40)

    dropout_rate = 0.15  # 15% expected dropout
    adjusted_n = math.ceil(ss_result.n_per_arm / (1 - dropout_rate))

    print(f"Expected dropout: {dropout_rate:.0%}")
    print(f"Unadjusted N:     {ss_result.n_per_arm}")
    print(f"Adjusted N:       {adjusted_n}")
    print(f"Total N:          {adjusted_n * 2}")

    # 4. Power table for different sample sizes
    print("\n" + "=" * 60)
    print("4. POWER TABLE (effect=0.5, alpha=0.05)")
    print("-" * 40)
    print("N per arm    Power")
    print("-" * 40)

    for n in [30, 40, 50, 60, 70, 80, 90, 100]:
        pwr = trial.estimate_power_analytical(
            n_per_arm=n,
            effect_size=0.5,
            sd=1.0,
            alpha=0.05
        )
        print(f"   {n:3d}        {pwr.power:.1%}")

    # 5. Sample size table for different effect sizes
    print("\n" + "=" * 60)
    print("5. SAMPLE SIZE TABLE (80% power, alpha=0.05)")
    print("-" * 40)
    print("Effect Size      N per arm    Total N")
    print("-" * 40)

    for effect, label in [(0.2, "small"), (0.3, ""), (0.4, ""),
                          (0.5, "medium"), (0.6, ""), (0.7, ""),
                          (0.8, "large")]:
        ss = trial.estimate_sample_size(
            target_power=0.80,
            effect_size=effect,
            alpha=0.05
        )
        label_str = f" ({label})" if label else ""
        print(f"   {effect:.1f}{label_str:10s}    {ss.n_per_arm:5d}        {ss.total_n:5d}")

    # 6. Effect size interpretation
    print("\n" + "=" * 60)
    print("6. EFFECT SIZE INTERPRETATION")
    print("-" * 40)
    print("Cohen's d | Overlap | Interpretation")
    print("-" * 40)
    print("   0.2    |  85%    | Small - subtle, often not visible")
    print("   0.5    |  67%    | Medium - noticeable difference")
    print("   0.8    |  53%    | Large - obvious difference")

    return {
        "power_result": result,
        "sample_size_result": ss_result,
        "adjusted_n": adjusted_n
    }


if __name__ == "__main__":
    main()
