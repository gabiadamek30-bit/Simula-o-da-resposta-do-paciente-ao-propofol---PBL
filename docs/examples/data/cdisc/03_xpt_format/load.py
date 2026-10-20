#!/usr/bin/env python3
"""
CDISC XPT Format Import - Python Example

Run: python load.py

Note: Requires pyreadstat or pandas with SAS support.
This example demonstrates the API - actual XPT files would be binary.
"""

from neopkpd.data import load_cdisc_xpt


def main():
    print("CDISC XPT Format Import")
    print("=" * 50)

    # Example usage (would need actual XPT files)
    print("\nAPI Usage:")
    print("-" * 50)
    print("""
    from neopkpd.data import load_cdisc_xpt

    # Load from XPT files (SAS transport format)
    data = load_cdisc_xpt(
        pc="pc.xpt",
        ex="ex.xpt",
        dm="dm.xpt",
        vs="vs.xpt"  # Optional vital signs for covariates
    )

    # Access loaded data
    print(f"Subjects: {len(data.subjects)}")
    print(f"Observations: {len(data.observations)}")

    # Convert to NeoPKPD format
    neopkpd_data = data.to_neopkpd()
    """)

    # XPT vs CSV comparison
    print("\nXPT vs CSV Format:")
    print("-" * 50)
    print("| Feature         | XPT              | CSV              |")
    print("|-----------------|------------------|------------------|")
    print("| Format          | Binary           | Text             |")
    print("| Size            | Smaller          | Larger           |")
    print("| Metadata        | Full SAS labels  | Headers only     |")
    print("| Regulatory use  | Standard         | Accepted         |")
    print("| Human readable  | No               | Yes              |")

    print("\nNote: XPT files are binary. This example shows the API.")
    print("For actual XPT files, use load_cdisc_xpt() function.")


if __name__ == "__main__":
    main()
