# Example 02: Oral First-Order Absorption
#
# This example demonstrates the OneCompOralFirstOrder model
# with first-order absorption kinetics.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "core", "NeoPKPD"))

using NeoPKPD

# Model parameters
# Ka = 1.5 /h (absorption rate constant)
# CL = 5.0 L/h (clearance)
# V = 50.0 L (volume of distribution)
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)

# 200 mg oral dose at t=0
doses = [DoseEvent(0.0, 200.0)]

# Create model specification
spec = ModelSpec(OneCompOralFirstOrder(), "oral_example", params, doses)

# Simulation from 0 to 24 hours with 0.25h resolution
grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))

# High-precision solver settings
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run simulation
result = simulate(spec, grid, solver)

# Calculate PK metrics
t = result.t
c = result.observations[:conc]

cmax_val = cmax(t, c)
auc_val = auc_trapezoid(t, c)

# Find Tmax (time of maximum concentration)
tmax_idx = argmax(c)
tmax_val = t[tmax_idx]

println("=== Oral First-Order Absorption Results ===")
println("Cmax: $(round(cmax_val, digits=3)) mg/L")
println("Tmax: $(round(tmax_val, digits=2)) h")
println("AUC(0-24): $(round(auc_val, digits=2)) mg·h/L")
println()
println("Time (h) | Concentration (mg/L)")
println("-" ^ 35)
for i in 1:4:length(t)
    println("$(lpad(t[i], 6)) | $(round(c[i], digits=4))")
end

# Write artifact
outpath = joinpath(@__DIR__, "..", "output", "02_oral_absorption.json")
write_execution_json(
    outpath;
    model_spec=spec,
    grid=grid,
    solver=solver,
    result=result
)
println("\nWrote: $outpath")
