#!/usr/bin/env python3
"""
Example: PK IV Bolus Simulation

This example demonstrates how to run a simple one-compartment
IV bolus PK simulation and compute exposure metrics.
"""

import neopkpd

# Initialize Julia (required once per session)
neopkpd.init_julia()

print("NeoPKPD version:", neopkpd.version())

# Run a simple IV bolus simulation
result = neopkpd.simulate_pk_iv_bolus(
    cl=1.0,
    v=10.0,
    doses=[{"time": 0.0, "amount": 100.0}],
    t0=0.0,
    t1=24.0,
    saveat=[0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0],
)

print("\nSimulation complete!")
print(f"Time points: {len(result['t'])}")
print(f"States: {list(result['states'].keys())}")
print(f"Observations: {list(result['observations'].keys())}")

# Compute metrics
print("\nPK Metrics:")
print(f"  Cmax: {neopkpd.cmax(result):.4f}")
print(f"  AUC: {neopkpd.auc_trapezoid(result):.4f}")
print(f"  Half-life: {neopkpd.half_life(1.0, 10.0):.4f} hours")

# Show concentration profile
print("\nConcentration profile:")
for t, c in zip(result['t'], result['observations']['conc']):
    print(f"  t={t:5.1f}h: {c:.4f}")

print("\nExample completed successfully!")
