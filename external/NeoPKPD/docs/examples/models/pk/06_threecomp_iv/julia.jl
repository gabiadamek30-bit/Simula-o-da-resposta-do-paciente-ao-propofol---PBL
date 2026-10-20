# Three-Compartment IV Bolus Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Three-Compartment IV Bolus Model")
println("="^50)

# Model parameters
CL = 5.0    # Clearance (L/h)
V1 = 10.0   # Central volume (L)
Q2 = 20.0   # Shallow distribution clearance (L/h)
V2 = 20.0   # Shallow peripheral volume (L)
Q3 = 2.0    # Deep distribution clearance (L/h)
V3 = 100.0  # Deep peripheral volume (L)
Dose = 100.0  # mg

# Create model specification
model = create_model_spec(
    "ThreeCompIVBolus";
    name = "threecomp_iv_example",
    params = Dict("CL" => CL, "V1" => V1, "Q2" => Q2, "V2" => V2, "Q3" => Q3, "V3" => V3),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings - need long time for deep compartment equilibration
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 168.0,  # 7 days
    saveat = vcat(collect(0.0:0.05:1.0), collect(1.0:0.5:24.0), collect(24.0:4.0:168.0))
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
auc = sum(0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i]) for i in 1:(length(times)-1))

# Derived parameters
Vss = V1 + V2 + V3

println("\nResults:")
println("-"^50)
println("Volumes:")
println("  V1 (central) = $(V1) L")
println("  V2 (shallow) = $(V2) L")
println("  V3 (deep)    = $(V3) L")
println("  Vss (total)  = $(Vss) L")
println("\nDistribution clearances:")
println("  Q2 (shallow) = $(Q2) L/h")
println("  Q3 (deep)    = $(Q3) L/h")
println("  CL (elim)    = $(CL) L/h")
println("\nSimulated PK:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")

# Show tri-exponential profile
println("\nConcentration-Time Profile (tri-exponential):")
println("-"^50)
println("Time (h)    Conc (mg/L)   Phase")
println("-"^50)
for (t, phase) in [(0.0, "peak"), (0.1, "α (rapid dist)"), (0.5, "α (rapid dist)"),
                   (1.0, "α→β transition"), (4.0, "β (slow dist)"),
                   (12.0, "β (slow dist)"), (24.0, "β→γ transition"),
                   (48.0, "γ (terminal)"), (96.0, "γ (terminal)"), (168.0, "γ (terminal)")]
    idx = argmin(abs.(times .- t))
    println("  $(lpad(string(t), 5))       $(lpad(string(round(conc[idx], digits=5)), 9))    $phase")
end
