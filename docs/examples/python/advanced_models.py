#!/usr/bin/env python3
"""
NeoPKPD Advanced Models Example

This example demonstrates the use of advanced PK and PD models available
in NeoPKPD, including:

PK Models:
- Two-compartment IV bolus (TwoCompIVBolus)
- Two-compartment oral (TwoCompOral)
- Three-compartment IV bolus (ThreeCompIVBolus)
- Transit compartment absorption (TransitAbsorption)
- Michaelis-Menten elimination (MichaelisMentenElimination)

PD Models:
- Sigmoid Emax (Hill equation)
- Biophase equilibration (effect compartment)
"""

import neopkpd


def run_twocomp_iv():
    """Two-compartment IV bolus simulation."""
    print("\n" + "=" * 60)
    print("Two-Compartment IV Bolus Model")
    print("=" * 60)

    result = neopkpd.simulate_pk_twocomp_iv_bolus(
        cl=10.0,      # Clearance (L/h)
        v1=50.0,      # Central volume (L)
        q=5.0,        # Inter-compartmental clearance (L/h)
        v2=100.0,     # Peripheral volume (L)
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0,
        t1=48.0,
        saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0],
    )

    print(f"Time points: {result['t'][:5]}...")
    print(f"Concentrations: {[f'{c:.2f}' for c in result['observations']['conc'][:5]]}...")
    print(f"Cmax: {neopkpd.cmax(result):.2f} mg/L")
    print(f"AUC (0-48h): {neopkpd.auc_trapezoid(result):.2f} mg*h/L")


def run_twocomp_oral():
    """Two-compartment oral absorption simulation."""
    print("\n" + "=" * 60)
    print("Two-Compartment Oral Model")
    print("=" * 60)

    result = neopkpd.simulate_pk_twocomp_oral(
        ka=1.5,       # Absorption rate constant (1/h)
        cl=10.0,      # Clearance (L/h)
        v1=50.0,      # Central volume (L)
        q=5.0,        # Inter-compartmental clearance (L/h)
        v2=100.0,     # Peripheral volume (L)
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0,
        t1=48.0,
        saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0],
    )

    print(f"Cmax: {neopkpd.cmax(result):.2f} mg/L")
    print(f"Tmax: {neopkpd.tmax(result):.2f} h")
    print(f"AUC (0-48h): {neopkpd.auc_trapezoid(result):.2f} mg*h/L")


def run_threecomp_iv():
    """Three-compartment IV bolus simulation."""
    print("\n" + "=" * 60)
    print("Three-Compartment IV Bolus Model")
    print("=" * 60)

    result = neopkpd.simulate_pk_threecomp_iv_bolus(
        cl=10.0,      # Clearance (L/h)
        v1=20.0,      # Central volume (L)
        q2=20.0,      # Shallow peripheral clearance (L/h)
        v2=50.0,      # Shallow peripheral volume (L)
        q3=2.0,       # Deep peripheral clearance (L/h)
        v3=200.0,     # Deep peripheral volume (L)
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0,
        t1=72.0,
        saveat=[0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0, 72.0],
    )

    print(f"Cmax: {neopkpd.cmax(result):.2f} mg/L")
    print(f"AUC (0-72h): {neopkpd.auc_trapezoid(result):.2f} mg*h/L")
    print("Note: Three-compartment shows tri-exponential decline")


def run_transit_absorption():
    """Transit compartment absorption simulation."""
    print("\n" + "=" * 60)
    print("Transit Compartment Absorption Model")
    print("=" * 60)

    result = neopkpd.simulate_pk_transit_absorption(
        n=5,          # Number of transit compartments
        ktr=1.5,      # Transit rate constant (1/h)
        ka=0.8,       # Absorption rate constant (1/h)
        cl=10.0,      # Clearance (L/h)
        v=100.0,      # Volume of distribution (L)
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0,
        t1=24.0,
        saveat=[0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0],
    )

    print(f"Cmax: {neopkpd.cmax(result):.2f} mg/L")
    print(f"Tmax: {neopkpd.tmax(result):.2f} h (delayed due to transit)")
    print(f"AUC (0-24h): {neopkpd.auc_trapezoid(result):.2f} mg*h/L")


def run_michaelis_menten():
    """Michaelis-Menten (saturable) elimination simulation."""
    print("\n" + "=" * 60)
    print("Michaelis-Menten Elimination Model")
    print("=" * 60)

    result = neopkpd.simulate_pk_michaelis_menten(
        vmax=15.0,    # Maximum elimination rate (mg/h)
        km=2.5,       # Michaelis constant (mg/L)
        v=50.0,       # Volume of distribution (L)
        doses=[{"time": 0.0, "amount": 500.0}],
        t0=0.0,
        t1=48.0,
        saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0],
    )

    print(f"Cmax: {neopkpd.cmax(result):.2f} mg/L")
    print(f"AUC (0-48h): {neopkpd.auc_trapezoid(result):.2f} mg*h/L")
    print("Note: At high concentrations, elimination is saturable (zero-order)")
    print("      At low concentrations, elimination becomes first-order")


def run_sigmoid_emax():
    """Sigmoid Emax (Hill equation) PD simulation."""
    print("\n" + "=" * 60)
    print("Sigmoid Emax (Hill Equation) PD Model")
    print("=" * 60)

    result = neopkpd.simulate_pkpd_sigmoid_emax(
        cl=5.0,       # PK clearance (L/h)
        v=50.0,       # PK volume (L)
        doses=[{"time": 0.0, "amount": 100.0}],
        e0=0.0,       # Baseline effect
        emax=100.0,   # Maximum effect
        ec50=1.0,     # Concentration for 50% effect (mg/L)
        gamma=2.0,    # Hill coefficient (sigmoidicity)
        t0=0.0,
        t1=24.0,
        saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
    )

    print(f"Peak concentration: {neopkpd.cmax(result):.2f} mg/L")
    print(f"Peak effect: {max(result['observations']['effect']):.1f}")
    print("Note: gamma > 1 creates steeper dose-response curve")


def run_biophase_equilibration():
    """Biophase equilibration (effect compartment) PD simulation."""
    print("\n" + "=" * 60)
    print("Biophase Equilibration (Effect Compartment) PD Model")
    print("=" * 60)

    result = neopkpd.simulate_pkpd_biophase_equilibration(
        cl=5.0,       # PK clearance (L/h)
        v=50.0,       # PK volume (L)
        doses=[{"time": 0.0, "amount": 100.0}],
        ke0=0.3,      # Effect compartment equilibration rate (1/h)
        e0=10.0,      # Baseline effect
        emax=90.0,    # Maximum effect
        ec50=1.0,     # Effect site conc. for 50% effect (mg/L)
        t0=0.0,
        t1=24.0,
        saveat=[0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0],
    )

    print(f"Peak concentration: {neopkpd.cmax(result):.2f} mg/L")
    print(f"Peak effect: {max(result['observations']['effect']):.1f}")
    print("Note: Effect lags behind concentration due to biophase equilibration")


def main():
    """Run all advanced model examples."""
    print("NeoPKPD Advanced Models Example")
    print("================================")

    # Initialize Julia once
    neopkpd.init_julia()

    # Run all examples
    run_twocomp_iv()
    run_twocomp_oral()
    run_threecomp_iv()
    run_transit_absorption()
    run_michaelis_menten()
    run_sigmoid_emax()
    run_biophase_equilibration()

    print("\n" + "=" * 60)
    print("All advanced model examples completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
