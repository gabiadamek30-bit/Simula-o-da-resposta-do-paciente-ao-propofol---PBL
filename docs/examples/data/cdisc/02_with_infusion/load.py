#!/usr/bin/env python3
"""
CDISC Infusion Data Import - Python Example

Run: python load.py
"""

from neopkpd.data import load_cdisc


def main():
    print("CDISC Infusion Data Import")
    print("=" * 50)

    # Load CDISC domains with infusion data
    data = load_cdisc(
        pc="../01_basic_pc_ex/pc.csv",  # Reuse PC data
        ex="ex.csv"  # EX with EXDUR
    )

    print(f"\nData Summary:")
    print(f"  Subjects: {len(data.subjects)}")
    print(f"  Doses:    {len(data.doses)}")

    # Display dosing with infusion duration
    print("\nInfusion Dosing:")
    print("-" * 50)
    print("Subject      Amount    Duration    Route")
    print("-" * 50)
    for dose in data.doses:
        duration_str = f"{dose.duration:.1f} h" if dose.duration else "Bolus"
        print(f"  {dose.subject_id}   {dose.amount:6.0f} mg   {duration_str:10s}   {dose.route}")

    # Convert to NeoPKPD format
    neopkpd_data = data.to_neopkpd()

    print("\nNeoPKPD Dose Events:")
    print("-" * 50)
    for subj in neopkpd_data.subjects[:1]:
        print(f"Subject {subj.id}:")
        for dose in subj.doses:
            if hasattr(dose, 'duration') and dose.duration:
                print(f"  Infusion: {dose.amount} mg over {dose.duration} h at t={dose.time}h")
            else:
                print(f"  Bolus: {dose.amount} mg at t={dose.time}h")

    return data


if __name__ == "__main__":
    main()
