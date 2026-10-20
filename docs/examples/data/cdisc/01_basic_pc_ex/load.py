#!/usr/bin/env python3
"""
Basic CDISC Data Import - Python Example

Run: python load.py
"""

from neopkpd.data import load_cdisc


def main():
    print("CDISC PC/EX Data Import")
    print("=" * 50)

    # Load CDISC domains
    data = load_cdisc(
        pc="pc.csv",
        ex="ex.csv",
        dm="dm.csv"
    )

    print(f"\nData Summary:")
    print(f"  Subjects:     {len(data.subjects)}")
    print(f"  Observations: {len(data.observations)}")
    print(f"  Doses:        {len(data.doses)}")

    # Display subject information
    print("\nSubjects:")
    print("-" * 50)
    for subj in data.subjects:
        print(f"  {subj.id}:")
        print(f"    Age: {subj.covariates.get('age', 'N/A')}")
        print(f"    Sex: {subj.covariates.get('sex', 'N/A')}")
        print(f"    Observations: {len([o for o in data.observations if o.subject_id == subj.id])}")

    # Display observations for first subject
    print("\nObservations (Subject 1):")
    print("-" * 50)
    print("Time (h)  Concentration (ng/mL)")
    print("-" * 50)
    subj1_obs = [o for o in data.observations if o.subject_id == data.subjects[0].id]
    for obs in subj1_obs:
        print(f"  {obs.time:5.1f}     {obs.dv:8.2f}")

    # Display dosing
    print("\nDosing:")
    print("-" * 50)
    for dose in data.doses:
        print(f"  Subject {dose.subject_id}: {dose.amount} mg at t={dose.time}h ({dose.route})")

    # Convert to NeoPKPD format for analysis
    print("\n" + "=" * 50)
    print("NeoPKPD Format")
    print("=" * 50)

    neopkpd_data = data.to_neopkpd()
    print(f"\nReady for analysis with {len(neopkpd_data.subjects)} subjects")

    return data


if __name__ == "__main__":
    main()
