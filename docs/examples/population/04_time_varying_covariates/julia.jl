# Time-Varying Covariates - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD
using Random

println("Time-Varying Covariates (Renal Function)")
println("="^50)

Random.seed!(12345)

# Example: Patient with declining renal function
# CrCL: 100 mL/min at t=0, declining to 50 mL/min by day 7

n_subjects = 50

# Generate time-varying CrCL for each subject
function generate_crcl_trajectory(baseline_crcl, decline_rate, times)
    return [max(baseline_crcl - decline_rate * t, 20.0) for t in times]
end

# Base model
base_model = create_model_spec(
    "OneCompIVBolus";
    name = "time_varying_cov_example",
    params = Dict("CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Time points where covariates are measured
covariate_times = collect(0.0:24.0:168.0)  # Daily for 7 days

# Generate covariate data
covariates_data = []
for i in 1:n_subjects
    baseline_crcl = 100.0 + 20.0 * randn()  # Variable baseline
    decline_rate = max(0, 7.0 + 3.0 * randn())  # Variable decline
    crcl_values = generate_crcl_trajectory(baseline_crcl, decline_rate, covariate_times)
    push!(covariates_data, Dict(
        :CRCL => Dict(zip(covariate_times, crcl_values))  # Time-varying
    ))
end

# Covariate model with time-varying effect
covariate_model = create_covariate_model(
    effects = [
        CovariateEffect(:CRCL, :CL, :power, exponent=0.5, reference=100.0, time_varying=true)
    ]
)

# Population specification
pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIV(omegas = Dict("CL" => 0.2, "V" => 0.15), seed = 12345),
    covariate_model = covariate_model,
    covariates = covariates_data,
    n = n_subjects
)

# Simulation settings
grid = SimulationGrid(t0 = 0.0, t1 = 168.0, saveat = collect(0.0:2.0:168.0))
solver = SolverConfig(alg = "Tsit5", reltol = 1e-8, abstol = 1e-10)

println("\nSimulating $(n_subjects) subjects over 7 days...")
println("CrCL declines over time, affecting drug clearance")

result = simulate_population(pop_spec, grid, solver)

println("\nExample Subject CrCL Trajectory:")
println("-"^40)
println("Day   CrCL (mL/min)  CL adjustment")
println("-"^40)
example_crcl = covariates_data[1][:CRCL]
for t in [0.0, 48.0, 96.0, 168.0]
    crcl = get(example_crcl, t, example_crcl[maximum(keys(example_crcl))])
    cl_adj = (crcl / 100.0)^0.5
    println("  $(Int(t/24))     $(round(crcl, digits=1))         $(round(cl_adj, digits=3))")
end

println("\nNote: Drug accumulation may increase as CrCL decreases")
