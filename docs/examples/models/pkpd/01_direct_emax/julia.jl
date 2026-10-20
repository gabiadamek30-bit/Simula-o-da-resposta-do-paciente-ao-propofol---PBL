# Direct Emax PD Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Direct Emax PD Model")
println("="^50)

# PK parameters
CL = 5.0   # L/h
V = 50.0   # L
Dose = 100.0  # mg

# PD parameters
E0 = 0.0      # Baseline effect
Emax = 100.0  # Maximum effect
EC50 = 1.0    # mg/L

# Create PKPD model specification
model = create_model_spec(
    "DirectEmax";
    name = "direct_emax_example",
    params = Dict(
        "CL" => CL, "V" => V,
        "E0" => E0, "Emax" => Emax, "EC50" => EC50
    ),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 48.0,
    saveat = collect(0.0:0.25:48.0)
)

solver = SolverConfig(
    alg = "Tsit5",
    reltol = 1e-10,
    abstol = 1e-12
)

# Run simulation
println("\nRunning PKPD simulation...")
result = simulate(model, grid, solver)

# Extract results
times = result.t
conc = result.observations["conc"]
effect = result.observations["effect"]

# Compute metrics
cmax = maximum(conc)
emax_observed = maximum(effect)
e_at_ec50 = E0 + Emax / 2  # Theoretical effect at EC50

println("\nResults:")
println("-"^50)
println("PK Parameters: CL=$(CL) L/h, V=$(V) L")
println("PD Parameters: E0=$(E0), Emax=$(Emax), EC50=$(EC50) mg/L")
println("\nSimulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Max Effect = $(round(emax_observed, digits=2))")
println("  Effect at C=EC50 (theoretical) = $(round(e_at_ec50, digits=2))")

# Show PK/PD relationship
println("\nConcentration-Effect Profile:")
println("-"^45)
println("Time (h)  Conc (mg/L)  Effect")
println("-"^45)
for t in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0]
    idx = argmin(abs.(times .- t))
    c = conc[idx]
    e = effect[idx]
    println("  $(lpad(string(t), 4))     $(lpad(string(round(c, digits=3)), 7))    $(round(e, digits=2))")
end

# Demonstrate no hysteresis (effect tracks concentration)
println("\nNote: Direct effect model - no hysteresis between PK and PD")
println("Effect at any time = E0 + Emax × C / (EC50 + C)")
