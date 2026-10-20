# Sigmoid Emax (Hill) PD Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Sigmoid Emax (Hill) PD Model")
println("="^50)

# PK parameters
CL = 5.0   # L/h
V = 50.0   # L
Dose = 100.0  # mg

# PD parameters
E0 = 0.0      # Baseline effect
Emax = 100.0  # Maximum effect
EC50 = 1.0    # mg/L
gamma = 2.0   # Hill coefficient

# Create PKPD model specification
model = create_model_spec(
    "SigmoidEmax";
    name = "sigmoid_emax_example",
    params = Dict(
        "CL" => CL, "V" => V,
        "E0" => E0, "Emax" => Emax, "EC50" => EC50, "gamma" => gamma
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

println("\nResults:")
println("-"^50)
println("PK Parameters: CL=$(CL) L/h, V=$(V) L")
println("PD Parameters:")
println("  E0 = $(E0), Emax = $(Emax), EC50 = $(EC50) mg/L")
println("  gamma (Hill coefficient) = $(gamma)")
println("\nSimulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Max Effect = $(round(emax_observed, digits=2))")

# Show steepness effect
println("\nEffect of Hill Coefficient (gamma = $(gamma)):")
println("  At C = 0.5 × EC50: Effect = $(round(E0 + Emax * (0.5^gamma) / (1 + 0.5^gamma), digits=1))%")
println("  At C = EC50:       Effect = $(round(E0 + Emax * 0.5, digits=1))%")
println("  At C = 2 × EC50:   Effect = $(round(E0 + Emax * (2^gamma) / (1 + 2^gamma), digits=1))%")

# Concentration-Effect profile
println("\nConcentration-Effect Profile:")
println("-"^45)
println("Time (h)  Conc (mg/L)  Effect")
println("-"^45)
for t in [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    idx = argmin(abs.(times .- t))
    println("  $(lpad(string(t), 4))     $(lpad(string(round(conc[idx], digits=3)), 7))    $(round(effect[idx], digits=2))")
end
