#!/usr/bin/env python3
"""
Monolix pk1cpt Import - Python Example

Run: python convert.py
"""

from neopkpd import import_monolix, simulate


def main():
    print("Monolix pk1cpt Import (One-compartment Oral)")
    print("=" * 50)

    # Import Monolix project
    model = import_monolix("project.mlxtran")

    print(f"\nImported Model:")
    print(f"  Type:       {model.model_type}")
    print(f"  Parameters: {model.parameters}")

    spec = model.spec
    print("\nNeoPKPD Specification:")
    print(f"  Model:   {spec.model}")
    print(f"  Ka:      {spec.params['Ka']} 1/h")
    print(f"  CL:      {spec.params['CL']} L/h")
    print(f"  V:       {spec.params['V']} L")

    if spec.iiv:
        print(f"  IIV omegas: {spec.iiv.omegas}")

    if spec.error:
        print(f"  Error:   {spec.error}")

    # Simulate
    result = simulate(
        spec,
        t_end=24.0,
        saveat=0.5,
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    print(f"\nSimulation:")
    print(f"  Time points: {len(result['times'])}")
    print(f"  Cmax: {max(result['observations']['conc']):.3f} mg/L")

    return model


if __name__ == "__main__":
    main()
