# Static Covariates - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD
using Random

println("Static Covariates (Weight, Age)")
println("="^50)

Random.seed!(12345)

# Generate covariate data for 100 subjects
n_subjects = 100
weights = clamp.(70.0 .+ 15.0 .* randn(n_subjects), 40.0, 120.0)
ages = clamp.(50.0 .+ 15.0 .* randn(n_subjects), 18.0, 80.0)

# Base model
base_model = create_model_spec(
    "OneCompIVBolus";
    name = "covariate_example",
    params = Dict("CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Covariate model
covariate_model = create_covariate_model(
    effects = [
        CovariateEffect(:WT, :CL, :power, exponent=0.75, reference=70.0),
        CovariateEffect(:WT, :V, :power, exponent=1.0, reference=70.0),
        CovariateEffect(:AGE, :CL, :linear, slope=-0.01, reference=40.0)
    ]
)

# Create covariates for each subject
covariates = [Dict(:WT => weights[i], :AGE => ages[i]) for i in 1:n_subjects]

# Population specification
pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIV(
        omegas = Dict("CL" => 0.25, "V" => 0.15),  # Reduced after covariate inclusion
        seed = 12345
    ),
    covariate_model = covariate_model,
    covariates = covariates,
    n = n_subjects
)

# Simulation settings
grid = SimulationGrid(t0 = 0.0, t1 = 24.0, saveat = collect(0.0:0.5:24.0))
solver = SolverConfig(alg = "Tsit5", reltol = 1e-8, abstol = 1e-10)

# Run simulation
println("\nRunning population simulation with covariates...")
result = simulate_population(pop_spec, grid, solver)

# Analyze covariate effects
println("\nCovariate Distribution:")
println("-"^40)
println("Weight: mean=$(round(mean(weights), digits=1)) kg, SD=$(round(std(weights), digits=1))")
println("Age:    mean=$(round(mean(ages), digits=1)) y, SD=$(round(std(ages), digits=1))")

# Show individual parameters
individual_params = result.params
CL_values = [p["CL"] for p in individual_params]
V_values = [p["V"] for p in individual_params]

println("\nIndividual Parameters (after covariate adjustment):")
println("-"^50)
println("CL: mean=$(round(mean(CL_values), digits=2)) L/h, range=[$(round(minimum(CL_values), digits=2)), $(round(maximum(CL_values), digits=2))]")
println("V:  mean=$(round(mean(V_values), digits=1)) L, range=[$(round(minimum(V_values), digits=1)), $(round(maximum(V_values), digits=1))]")

# Example calculations
println("\nExample Covariate Effects:")
println("-"^50)
println("50 kg, 30y: CL = 5.0 × (50/70)^0.75 × (1-0.01×(30-40)) = $(round(5.0 * (50/70)^0.75 * (1 - 0.01*(30-40)), digits=2)) L/h")
println("90 kg, 70y: CL = 5.0 × (90/70)^0.75 × (1-0.01×(70-40)) = $(round(5.0 * (90/70)^0.75 * (1 - 0.01*(70-40)), digits=2)) L/h")
