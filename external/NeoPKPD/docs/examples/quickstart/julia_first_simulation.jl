# NeoPKPD Quickstart - Julia
# Run: julia --project=packages/core docs/examples/quickstart/julia_first_simulation.jl

using NeoPKPD

println("NeoPKPD Quickstart - Julia")
println("="^40)

# 1. Create a one-compartment IV bolus model
model = create_model_spec(
    "OneCompIVBolus";
    name = "quickstart_example",
    params = Dict("CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

println("\nModel: $(model.kind)")
println("Parameters: CL = $(model.params["CL"]) L/h, V = $(model.params["V"]) L")
println("Dose: $(model.doses[1].amount) mg IV at t=0")

# 2. Configure simulation
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 24.0,
    saveat = collect(0.0:0.5:24.0)
)

solver = SolverConfig(
    alg = "Tsit5",
    reltol = 1e-8,
    abstol = 1e-10
)

# 3. Run simulation
println("\nRunning simulation...")
result = simulate(model, grid, solver)

# 4. Extract results
times = result.t
conc = result.observations["conc"]

# 5. Compute metrics
cmax = maximum(conc)
tmax_idx = argmax(conc)
tmax = times[tmax_idx]

# AUC by trapezoidal rule
auc = 0.0
for i in 1:(length(times)-1)
    auc += 0.5 * (conc[i] + conc[i+1]) * (times[i+1] - times[i])
end

# Terminal half-life (from last few points)
log_conc = log.(conc[end-4:end])
t_segment = times[end-4:end]
slope = (log_conc[end] - log_conc[1]) / (t_segment[end] - t_segment[1])
t_half = -log(2) / slope

# 6. Display results
println("\n" * "="^40)
println("RESULTS")
println("="^40)
println("Cmax:     $(round(cmax, digits=3)) mg/L")
println("Tmax:     $(round(tmax, digits=1)) h")
println("AUC₀₋₂₄:  $(round(auc, digits=2)) mg·h/L")
println("t½:       $(round(t_half, digits=2)) h")

# 7. Sample concentration-time data
println("\nConcentration-Time Profile (selected points):")
println("-"^30)
println("Time (h)    Conc (mg/L)")
println("-"^30)
for t in [0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]
    idx = findfirst(x -> x ≈ t, times)
    if !isnothing(idx)
        println("  $(lpad(string(t), 4))        $(round(conc[idx], digits=4))")
    end
end

println("\nQuickstart complete!")
