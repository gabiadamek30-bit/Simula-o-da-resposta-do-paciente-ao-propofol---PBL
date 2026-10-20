# Basic Inter-Individual Variability (IIV) - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Basic Inter-Individual Variability (IIV)")
println("="^50)

# Base PK model
base_model = create_model_spec(
    "OneCompIVBolus";
    name = "population_iiv_example",
    params = Dict("CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Population specification with IIV
n_subjects = 100
pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIV(
        omegas = Dict("CL" => 0.3, "V" => 0.2),  # Standard deviations
        seed = 12345  # For reproducibility
    ),
    n = n_subjects
)

# Simulation settings
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

# Run population simulation
println("\nRunning population simulation (n=$(n_subjects))...")
result = simulate_population(pop_spec, grid, solver)

# Extract individual parameters
individual_params = result.params
CL_values = [p["CL"] for p in individual_params]
V_values = [p["V"] for p in individual_params]

println("\nIndividual Parameter Distribution:")
println("-"^50)
println("Parameter  Mean    Median   SD      5th%    95th%")
println("-"^50)
println("CL (L/h)   $(round(mean(CL_values), digits=2))   $(round(median(CL_values), digits=2))    $(round(std(CL_values), digits=2))   $(round(quantile(CL_values, 0.05), digits=2))   $(round(quantile(CL_values, 0.95), digits=2))")
println("V (L)      $(round(mean(V_values), digits=1))   $(round(median(V_values), digits=1))   $(round(std(V_values), digits=1))   $(round(quantile(V_values, 0.05), digits=1))   $(round(quantile(V_values, 0.95), digits=1))")

# Population summary of concentrations
summaries = result.summaries
conc_summary = summaries["conc"]

println("\nPopulation Concentration Summary:")
println("-"^50)
println("Time (h)  Mean    Median  5th%    95th%")
println("-"^50)
times = result.t
for (i, t) in enumerate([0.0, 2.0, 4.0, 8.0, 12.0, 24.0])
    idx = findfirst(x -> x ≈ t, times)
    if !isnothing(idx)
        println("  $(lpad(string(t), 4))     $(round(conc_summary.mean[idx], digits=3))  $(round(conc_summary.median[idx], digits=3))  $(round(conc_summary.quantiles["0.05"][idx], digits=3))  $(round(conc_summary.quantiles["0.95"][idx], digits=3))")
    end
end

println("\nNote: IIV creates a distribution of PK profiles")
println("Each subject has unique CL and V values")
