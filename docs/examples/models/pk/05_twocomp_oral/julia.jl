# Two-Compartment Oral Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Two-Compartment Oral Model")
println("="^50)

# Model parameters
Ka = 1.5   # Absorption rate constant (1/h)
CL = 5.0   # Clearance (L/h)
V1 = 10.0  # Central volume (L)
Q = 10.0   # Inter-compartmental clearance (L/h)
V2 = 40.0  # Peripheral volume (L)
Dose = 100.0  # mg

# Create model specification
model = create_model_spec(
    "TwoCompOral";
    name = "twocomp_oral_example",
    params = Dict("Ka" => Ka, "CL" => CL, "V1" => V1, "Q" => Q, "V2" => V2),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 72.0,
    saveat = collect(0.0:0.25:72.0)
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

# Derived parameters
Vss = V1 + V2
k10 = CL / V1
k12 = Q / V1
k21 = Q / V2

println("\nResults:")
println("-"^50)
println("Parameters:")
println("  Ka = $(Ka) 1/h, CL = $(CL) L/h")
println("  V1 = $(V1) L, V2 = $(V2) L, Q = $(Q) L/h")
println("  Vss = $(Vss) L")
println("\nSimulated PK:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Tmax = $(round(tmax, digits=2)) h")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")

# Show triphasic profile
println("\nConcentration-Time Profile:")
println("-"^45)
for t in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0, 48.0, 72.0]
    idx = argmin(abs.(times .- t))
    println("t=$(lpad(string(t), 5))h: $(round(conc[idx], digits=4)) mg/L")
end
