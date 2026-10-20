#!/usr/bin/env python3
"""
NONMEM ADVAN1 Import - Python Example

Run: python convert.py
"""

from neopkpd import import_nonmem, simulate
import json


def main():
    print("NONMEM ADVAN1 Import")
    print("=" * 50)

    # Import NONMEM control file
    model = import_nonmem("run001.ctl")

    print(f"\nImported Model:")
    print(f"  Type:       {model.model_type}")
    print(f"  Parameters: {model.parameters}")

    # Display NeoPKPD specification
    print("\nNeoPKPD Specification:")
    print("-" * 50)
    spec = model.spec

    print(f"  Model:   {spec.model}")
    print(f"  Params:  {spec.params}")

    if spec.iiv:
        print(f"  IIV:")
        print(f"    Kind:   {spec.iiv.kind}")
        print(f"    Omegas: {spec.iiv.omegas}")

    if spec.error:
        print(f"  Error:   {spec.error}")

    # Compare with expected
    print("\n" + "=" * 50)
    print("VALIDATION")
    print("=" * 50)

    with open("expected.json") as f:
        expected = json.load(f)

    # Validate parameters
    params_match = spec.params == expected["params"]
    print(f"Parameters match: {params_match}")

    # Validate IIV
    iiv_match = spec.iiv.omegas == expected["iiv"]["omegas"]
    print(f"IIV match: {iiv_match}")

    # Simulate to verify
    print("\n" + "=" * 50)
    print("SIMULATION TEST")
    print("=" * 50)

    result = simulate(
        spec,
        t_end=24.0,
        saveat=1.0,
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    print(f"\nSimulation successful: {len(result['times'])} time points")
    print(f"Cmax: {max(result['observations']['conc']):.3f} mg/L")

    return model


if __name__ == "__main__":
    main()
