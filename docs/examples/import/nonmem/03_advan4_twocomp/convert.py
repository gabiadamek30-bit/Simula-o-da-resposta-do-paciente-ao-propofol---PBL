#!/usr/bin/env python3
"""
NONMEM ADVAN4 Import - Python Example

Run: python convert.py
"""

from neopkpd import import_nonmem, simulate


def main():
    print("NONMEM ADVAN4 Import (Two-compartment)")
    print("=" * 50)

    model = import_nonmem("run003.ctl")

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
    print(f"  Cmax (central):     {max(result['observations']['conc']):.3f} mg/L")
    print(f"  C at 24h:           {result['observations']['conc'][-1]:.4f} mg/L")

    return model


if __name__ == "__main__":
    main()
