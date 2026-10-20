# Complex Covariate Model - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD
using Random
using LinearAlgebra

println("Complex Covariate Model")
println("="^50)

Random.seed!(12345)

n_subjects = 200

# Generate covariates
weights = clamp.(70.0 .+ 15.0 .* randn(n_subjects), 40.0, 120.0)
ages = clamp.(50.0 .+ 15.0 .* randn(n_subjects), 18.0, 80.0)
crcl = clamp.(100.0 .+ 25.0 .* randn(n_subjects), 30.0, 150.0)
sex = rand(0:1, n_subjects)  # 0=male, 1=female

# Base model
base_model = create_model_spec(
    "OneCompOralFirstOrder";
    name = "complex_covariate_example",
    params = Dict("Ka" => 1.5, "CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Complex covariate model
covariate_model = create_covariate_model(
    effects = [
        CovariateEffect(:WT, :CL, :power, exponent=0.75, reference=70.0),
        CovariateEffect(:WT, :V, :power, exponent=1.0, reference=70.0),
        CovariateEffect(:CRCL, :CL, :power, exponent=0.5, reference=100.0),
        CovariateEffect(:SEX, :CL, :categorical, coefficient=0.85)  # Female = 85% of male
    ]
)

# Create covariates
covariates = [Dict(:WT => weights[i], :AGE => ages[i], :CRCL => crcl[i], :SEX => sex[i])
              for i in 1:n_subjects]

# Omega matrix with correlation between CL and V
omega_matrix = [0.09 0.02 0.0;
                0.02 0.04 0.0;
                0.0  0.0  0.16]

# Population specification with correlated random effects
pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIVCorrelated(
        parameters = ["CL", "V", "Ka"],
        omega_matrix = omega_matrix,
        seed = 12345
    ),
    covariate_model = covariate_model,
    covariates = covariates,
    n = n_subjects
)

# Simulation
grid = SimulationGrid(t0 = 0.0, t1 = 24.0, saveat = collect(0.0:0.5:24.0))
solver = SolverConfig(alg = "Tsit5", reltol = 1e-8, abstol = 1e-10)

println("\nRunning complex population simulation...")
println("  Subjects: $(n_subjects)")
println("  Covariates: WT, AGE, CrCL, SEX")
println("  Correlated IIV: CL-V correlation")

result = simulate_population(pop_spec, grid, solver)

# Analyze results
individual_params = result.params
CL_values = [p["CL"] for p in individual_params]
V_values = [p["V"] for p in individual_params]

println("\nPopulation Summary:")
println("-"^50)
println("Parameter   Mean    SD      CV%")
println("-"^50)
println("CL (L/h)    $(round(mean(CL_values), digits=2))   $(round(std(CL_values), digits=2))   $(round(std(CL_values)/mean(CL_values)*100, digits=1))")
println("V (L)       $(round(mean(V_values), digits=1))   $(round(std(V_values), digits=1))   $(round(std(V_values)/mean(V_values)*100, digits=1))")

# CL by sex
CL_male = [CL_values[i] for i in 1:n_subjects if sex[i] == 0]
CL_female = [CL_values[i] for i in 1:n_subjects if sex[i] == 1]
println("\nCL by Sex:")
println("  Male:   $(round(mean(CL_male), digits=2)) L/h (n=$(length(CL_male)))")
println("  Female: $(round(mean(CL_female), digits=2)) L/h (n=$(length(CL_female)))")
println("  Ratio:  $(round(mean(CL_female)/mean(CL_male), digits=2))")
