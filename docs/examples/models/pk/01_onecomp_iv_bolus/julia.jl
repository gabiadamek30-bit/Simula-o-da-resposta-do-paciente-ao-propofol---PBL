# One-Compartment IV Bolus Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("One-Compartment IV Bolus Model")
println("="^50)

# Model parameters
CL = 5.0   # Clearance (L/h)
V = 50.0   # Volume of distribution (L)
Dose = 100.0  # mg

# Create model specification
model = create_model_spec(
    "OneCompIVBolus";
    name = "onecomp_iv_bolus_example",
    params = Dict("CL" => CL, "V" => V),
    doses = [DoseEvent(time=0.0, amount=Dose)]
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 48.0,
    saveat = collect(0.0:0.5:48.0)
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
t_half_theoretical = log(2) / k
cmax_theoretical = Dose / V
auc_theoretical = Dose / CL

println("\nResults:")
println("-"^50)
println("Simulated:")
println("  Cmax = $(round(cmax, digits=4)) mg/L")
println("  Tmax = $(round(tmax, digits=2)) h")
println("  AUC  = $(round(auc, digits=2)) mg·h/L")
println("\nTheoretical:")
println("  Cmax = $(round(cmax_theoretical, digits=4)) mg/L")
println("  t½   = $(round(t_half_theoretical, digits=2)) h")
println("  AUC  = $(round(auc_theoretical, digits=2)) mg·h/L")

# Sample output
println("\nConcentration-Time Profile (first 12h):")
println("-"^30)
for t in 0.0:2.0:12.0
    idx = findfirst(x -> x ≈ t, times)
    if !isnothing(idx)
        println("t=$(lpad(string(Int(t)), 2))h: $(round(conc[idx], digits=4)) mg/L")
    end
end
