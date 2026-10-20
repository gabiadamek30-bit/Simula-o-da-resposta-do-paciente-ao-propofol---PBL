#!/usr/bin/env python3
"""
Monolix pk2cpt Import - Python Example

Run: python convert.py
"""

from neopkpd import import_monolix, simulate


def main():
    print("Monolix pk2cpt Import (Two-compartment IV)")
    print("=" * 50)

    model = import_monolix("project.mlxtran")

    print(f"\nImported Model:")
    print(f"  Type:       {model.model_type}")

    spec = model.spec
    print("\nNeoPKPD Specification:")
    print(f"  Model:   {spec.model}")
    print(f"  CL:      {spec.params['CL']} L/h")
    print(f"  V1:      {spec.params['V1']} L")
    print(f"  Q:       {spec.params['Q']} L/h")
    print(f"  V2:      {spec.params['V2']} L")

    # Simulate
    result = simulate(
        spec,
        t_end=24.0,
        saveat=0.1,
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    print(f"\nSimulation:")
    print(f"  Cmax: {max(result['observations']['conc']):.3f} mg/L")

    return model


if __name__ == "__main__":
    main()
