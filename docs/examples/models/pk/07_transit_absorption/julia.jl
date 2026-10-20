# Transit Compartment Absorption Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Transit Compartment Absorption Model")
println("="^50)

# Model parameters
Ktr = 2.0   # Transit rate constant (1/h)
n = 3       # Number of transit compartments
CL = 5.0    # Clearance (L/h)
V = 50.0    # Volume of distribution (L)
Dose = 100.0  # mg

# Calculate mean transit time
MTT = (n + 1) / Ktr

# Create model specification
model = create_model_spec(
    "TransitAbsorption";
    name = "transit_absorption_example",
    params = Dict("Ktr" => Ktr, "n" => n, "CL" => CL, "V" => V),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 48.0,
    saveat = collect(0.0:0.1:48.0)
)

solver = SolverConfig(
    alg = "Tsit5",
    reltol = 1e-10,
    abstol = 1e-12
)

# Run simulation
println("\nRunning simulation...")
result = simulate(model, grid, solver)

# Extract results
times = result.t
conc = result.observations["conc"]

# Compute metrics
cmax = maximum(conc)
tmax = times[argmax(conc)]
auc = sum(0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i]) for i in 1:(length(times)-1))

# Compare with first-order absorption
# For first-order with same MTT, Ka ≈ 1/MTT
Ka_equivalent = 1.0 / MTT

println("\nResults:")
println("-"^50)
println("Transit Model Parameters:")
println("  Ktr = $(Ktr) 1/h")
println("  n   = $(n) compartments")
println("  MTT = $(round(MTT, digits=2)) h")
println("\nPK Parameters:")
println("  CL = $(CL) L/h")
println("  V  = $(V) L")
println("\nSimulated PK:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Tmax = $(round(tmax, digits=2)) h")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")
println("\nNote: First-order model with equivalent MTT would have Ka ≈ $(round(Ka_equivalent, digits=2)) 1/h")

# Show absorption profile
println("\nConcentration-Time Profile:")
println("-"^40)
for t in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
    idx = argmin(abs.(times .- t))
    marker = t ≈ tmax ? " <-- Tmax" : ""
    println("t=$(lpad(string(t), 4))h: $(round(conc[idx], digits=4)) mg/L$marker")
end
