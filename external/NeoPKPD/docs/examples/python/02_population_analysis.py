#!/usr/bin/env python3
"""
Example 02: Population Analysis with Python

This example demonstrates population PK simulation and analysis
using the NeoPKPD Python bindings.
"""

import neopkpd

# Initialize Julia (required once per session)
print("Initializing Julia...")
neopkpd.init_julia()

# Check version
print(f"NeoPKPD version: {neopkpd.version()}")
print()

# Run population simulation
print("Running population simulation with 100 individuals...")
result = neopkpd.simulate_population_iv_bolus(
    cl=5.0,              # Typical clearance (L/h)
    v=50.0,              # Typical volume (L)
    doses=[{"time": 0.0, "amount": 100.0}],  # 100 mg IV bolus
    t0=0.0,
    t1=24.0,
    saveat=[float(t) * 0.5 for t in range(49)],  # Every 0.5h
    n=100,               # 100 individuals
    seed=12345,          # Reproducible
    omegas={
        "CL": 0.3,       # 30% CV on CL
        "V": 0.2         # 20% CV on V
    }
)

print(f"Simulation complete!")
print(f"Number of individuals: {len(result['individuals'])}")
print()

# Analyze parameter distributions
cls = [p["CL"] for p in result["params"]]
vs = [p["V"] for p in result["params"]]

print("=== Parameter Distributions ===")
print(f"CL: mean={sum(cls)/len(cls):.2f}, min={min(cls):.2f}, max={max(cls):.2f} L/h")
print(f"V:  mean={sum(vs)/len(vs):.2f}, min={min(vs):.2f}, max={max(vs):.2f} L")
print()

# Concentration summary
summary = result["summaries"]["conc"]
t = result["individuals"][0]["t"]

print("=== Concentration Summary (mg/L) ===")
print(f"{'Time (h)':<10} {'Mean':<10} {'Median':<10} {'5%':<10} {'95%':<10}")
print("-" * 50)

# Print at selected time points
for i in [0, 4, 8, 12, 24, 36, 48]:
    if i < len(t):
        print(f"{t[i]:<10.1f} {summary['mean'][i]:<10.3f} {summary['median'][i]:<10.3f} "
              f"{summary['quantiles']['0.05'][i]:<10.3f} {summary['quantiles']['0.95'][i]:<10.3f}")

print()

# Calculate individual Cmax values
cmax_values = [max(ind["observations"]["conc"]) for ind in result["individuals"]]
print("=== Individual Cmax Distribution ===")
print(f"Mean Cmax: {sum(cmax_values)/len(cmax_values):.3f} mg/L")
print(f"Min Cmax:  {min(cmax_values):.3f} mg/L")
print(f"Max Cmax:  {max(cmax_values):.3f} mg/L")

# Optional: Export to CSV (if pandas available)
try:
    import pandas as pd

    # Create long-format DataFrame
    data = []
    for i, ind in enumerate(result["individuals"]):
        for j, t_val in enumerate(ind["t"]):
            data.append({
                "id": i + 1,
                "time": t_val,
                "conc": ind["observations"]["conc"][j],
                "CL": result["params"][i]["CL"],
                "V": result["params"][i]["V"]
            })

    df = pd.DataFrame(data)
    print()
    print("=== DataFrame Summary ===")
    print(df.groupby("time")["conc"].describe().head(10))

except ImportError:
    print("\nNote: Install pandas for DataFrame export functionality")
