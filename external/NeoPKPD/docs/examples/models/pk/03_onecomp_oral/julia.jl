# One-Compartment Oral First-Order Absorption - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("One-Compartment Oral First-Order Absorption Model")
println("="^50)

# Model parameters
Ka = 1.5   # Absorption rate constant (1/h)
CL = 5.0   # Clearance (L/h)
V = 50.0   # Volume of distribution (L)
Dose = 100.0  # mg

# Create model specification
model = create_model_spec(
    "OneCompOralFirstOrder";
    name = "onecomp_oral_example",
    params = Dict("Ka" => Ka, "CL" => CL, "V" => V),
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
println("\nRunning simulation...")
result = simulate(model, grid, solver)

# Extract results
times = result.t
conc = result.observations["conc"]

# Compute metrics
cmax = maximum(conc)
tmax = times[argmax(conc)]
auc = sum(0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i]) for i in 1:(length(times)-1))

# Theoretical values
k = CL / V
tmax_theoretical = log(Ka/k) / (Ka - k)
cmax_theoretical = (Dose/V) * (Ka/(Ka-k)) * (exp(-k*tmax_theoretical) - exp(-Ka*tmax_theoretical))

println("\nResults:")
println("-"^50)
println("Parameters: Ka=$(Ka) 1/h, CL=$(CL) L/h, V=$(V) L")
println("\nSimulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Tmax = $(round(tmax, digits=2)) h")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")
println("\nTheoretical:")
println("  Cmax = $(round(cmax_theoretical, digits=4)) mg/L")
println("  Tmax = $(round(tmax_theoretical, digits=2)) h")

# Absorption vs elimination phases
println("\nConcentration-Time Profile:")
println("-"^40)
println("Time (h)  Conc (mg/L)  Phase")
println("-"^40)
for t in [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0]
    idx = argmin(abs.(times .- t))
    phase = t < tmax ? "absorption" : (t ≈ tmax ? "peak" : "elimination")
    println("  $(lpad(string(t), 4))     $(lpad(string(round(conc[idx], digits=4)), 8))    $phase")
end
