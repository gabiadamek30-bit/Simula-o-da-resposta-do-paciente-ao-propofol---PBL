#!/usr/bin/env python3
"""
NONMEM ADVAN2 Import - Python Example

Run: python convert.py
"""

from neopkpd import import_nonmem, simulate


def main():
    print("NONMEM ADVAN2 Import (Oral)")
    print("=" * 50)

    # Import NONMEM control file
    model = import_nonmem("run002.ctl")

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
        print(f"  IIV:")
        for param, omega in spec.iiv.omegas.items():
            cv = (omega ** 0.5) * 100
            print(f"    {param}: omega={omega:.4f} (CV~{cv:.0f}%)")

    # Simulate
    result = simulate(
        spec,
        t_end=24.0,
        saveat=0.5,
        doses=[{"time": 0.0, "amount": 100.0}]
    )

    print(f"\nSimulation:")
    print(f"  Cmax:  {max(result['observations']['conc']):.3f} mg/L")
    print(f"  Tmax:  {result['times'][result['observations']['conc'].index(max(result['observations']['conc']))]} h")

    return model


if __name__ == "__main__":
    main()
