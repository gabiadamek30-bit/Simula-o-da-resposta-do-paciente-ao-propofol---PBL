# Two-Compartment IV Bolus Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Two-Compartment IV Bolus Model")
println("="^50)

# Model parameters
CL = 5.0   # Clearance (L/h)
V1 = 10.0  # Central volume (L)
Q = 10.0   # Inter-compartmental clearance (L/h)
V2 = 40.0  # Peripheral volume (L)
Dose = 100.0  # mg

# Create model specification
model = create_model_spec(
    "TwoCompIVBolus";
    name = "twocomp_iv_example",
    params = Dict("CL" => CL, "V1" => V1, "Q" => Q, "V2" => V2),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 48.0,
    saveat = vcat(collect(0.0:0.05:1.0), collect(1.0:0.25:48.0))  # Dense early sampling
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

# Macro rate constants (eigenvalues)
sum_k = k10 + k12 + k21
diff = sqrt((k10 + k12 + k21)^2 - 4*k10*k21)
alpha = (sum_k + diff) / 2
beta = (sum_k - diff) / 2

t_half_alpha = log(2) / alpha
t_half_beta = log(2) / beta

println("\nResults:")
println("-"^50)
println("Micro-constants:")
println("  k10 = $(round(k10, digits=4)) 1/h")
println("  k12 = $(round(k12, digits=4)) 1/h")
println("  k21 = $(round(k21, digits=4)) 1/h")
println("\nMacro-constants:")
println("  α = $(round(alpha, digits=4)) 1/h (t½α = $(round(t_half_alpha, digits=2)) h)")
println("  β = $(round(beta, digits=4)) 1/h (t½β = $(round(t_half_beta, digits=2)) h)")
println("\nVolumes:")
println("  V1 = $(V1) L (central)")
println("  V2 = $(V2) L (peripheral)")
println("  Vss = $(Vss) L (steady-state)")
println("\nSimulated PK:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")

# Show bi-exponential profile
println("\nConcentration-Time Profile (bi-exponential):")
println("-"^45)
println("Time (h)  Conc (mg/L)   Phase")
println("-"^45)
for t in [0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0, 48.0]
    idx = argmin(abs.(times .- t))
    phase = t < 1.0 ? "distribution (α)" : "terminal (β)"
    println("  $(lpad(string(t), 5))     $(lpad(string(round(conc[idx], digits=4)), 8))    $phase")
end
