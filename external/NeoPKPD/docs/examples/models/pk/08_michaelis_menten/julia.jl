# Michaelis-Menten Elimination Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Michaelis-Menten Elimination Model")
println("="^50)

# Model parameters
Vmax = 50.0  # Maximum elimination rate (mg/h)
Km = 10.0    # Michaelis constant (mg/L)
V = 50.0     # Volume of distribution (L)
Dose = 500.0 # mg (high dose to show saturation)

# Calculate derived parameters
CLint = Vmax / Km  # Intrinsic clearance at low concentrations
t_half_linear = 0.693 * V * Km / Vmax

# Create model specification
model = create_model_spec(
    "MichaelisMentenElimination";
    name = "michaelis_menten_example",
    params = Dict("Vmax" => Vmax, "Km" => Km, "V" => V),
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
auc = sum(0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i]) for i in 1:(length(times)-1))

# Analyze saturation
C0 = Dose / V
println("\nResults:")
println("-"^50)
println("Model Parameters:")
println("  Vmax = $(Vmax) mg/h")
println("  Km   = $(Km) mg/L")
println("  V    = $(V) L")
println("\nDerived Parameters:")
println("  CLint (at low C) = $(round(CLint, digits=2)) L/h")
println("  t½ (at low C)    = $(round(t_half_linear, digits=2)) h")
println("\nSaturation Analysis:")
println("  C0 = $(round(C0, digits=2)) mg/L (Dose/V)")
println("  C0/Km = $(round(C0/Km, digits=1)) (>>1 = saturated)")
println("\nSimulated PK:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")

# Show elimination kinetics changing
println("\nConcentration-Time Profile (nonlinear elimination):")
println("-"^55)
println("Time (h)  Conc (mg/L)  C/Km    Kinetics")
println("-"^55)
for t in [0.0, 2.0, 4.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0]
    idx = argmin(abs.(times .- t))
    c = conc[idx]
    c_km = c / Km
    kinetics = c_km > 5 ? "zero-order (saturated)" : (c_km > 0.2 ? "mixed" : "first-order")
    println("  $(lpad(string(t), 4))      $(lpad(string(round(c, digits=3)), 7))    $(lpad(string(round(c_km, digits=2)), 5))   $kinetics")
end

# Compare to linear model
println("\nNote: A linear model with CL = CLint = $(round(CLint, digits=1)) L/h")
println("would have t½ = $(round(t_half_linear, digits=1)) h and AUC = $(round(Dose/CLint, digits=1)) mg·h/L")
println("The actual AUC is $(round(auc/(Dose/CLint) * 100 - 100, digits=0))% higher due to saturation.")
