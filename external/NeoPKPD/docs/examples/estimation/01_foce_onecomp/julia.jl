# FOCE-I One-Compartment Estimation - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD
using CSV
using DataFrames

println("FOCE-I One-Compartment Estimation")
println("="^50)

# Load observed data
data_path = joinpath(@__DIR__, "data.csv")
df = CSV.read(data_path, DataFrame)

println("\nData Summary:")
println("  Subjects: $(length(unique(df.ID)))")
println("  Observations: $(nrow(df))")
println("  Time range: $(minimum(df.TIME)) - $(maximum(df.TIME)) h")

# Create observed data object
observed = create_observed_data(
    df;
    id_col = :ID,
    time_col = :TIME,
    dv_col = :DV,
    amt_col = :AMT,
    evid_col = :EVID
)

# Model specification
model_spec = create_model_spec(
    "OneCompOralFirstOrder";
    name = "foce_example"
)

# Estimation configuration
config = EstimationConfig(
    method = FOCEI(),
    # Fixed effects (theta)
    theta_names = ["CL", "V", "Ka"],
    theta_init = [4.0, 40.0, 1.0],  # Initial guesses
    theta_lower = [0.1, 1.0, 0.1],
    theta_upper = [50.0, 500.0, 10.0],
    # Random effects (omega)
    omega_init = Dict("CL" => 0.1, "V" => 0.1, "Ka" => 0.2),
    omega_structure = :diagonal,
    # Residual error (sigma)
    error_model = :proportional,
    sigma_init = [0.05],
    # Algorithm settings
    maxiter = 500,
    abstol = 1e-6,
    reltol = 1e-4
)

# Run estimation
println("\nRunning FOCE-I estimation...")
result = estimate(observed, model_spec, config)

# Display results
println("\n" * "="^60)
println("ESTIMATION RESULTS")
println("="^60)

println("\nFixed Effects (θ):")
println("-"^50)
println("Parameter   Estimate    SE        RSE%")
println("-"^50)
for (name, est, se) in zip(config.theta_names, result.theta, result.theta_se)
    rse = se / est * 100
    println("$(rpad(name, 12))$(lpad(string(round(est, digits=3)), 8))    $(lpad(string(round(se, digits=3)), 6))    $(round(rse, digits=1))%")
end

println("\nRandom Effects (ω):")
println("-"^50)
for (name, omega) in result.omega
    println("ω_$(name) = $(round(sqrt(omega), digits=3)) ($(round(sqrt(omega)*100, digits=1))% CV)")
end

println("\nResidual Error (σ):")
println("-"^50)
println("σ_prop = $(round(result.sigma[1], digits=4)) ($(round(result.sigma[1]*100, digits=1))% CV)")

println("\nObjective Function:")
println("-"^50)
println("OFV = $(round(result.ofv, digits=2))")
println("AIC = $(round(result.aic, digits=2))")
println("BIC = $(round(result.bic, digits=2))")

println("\nConvergence:")
println("-"^50)
println("Status: $(result.converged ? "Converged" : "Not converged")")
println("Iterations: $(result.iterations)")

# Compare to true values
println("\n" * "="^60)
println("COMPARISON TO TRUE VALUES")
println("="^60)
println("Parameter   Estimate   True    Bias%")
println("-"^50)
true_values = [5.0, 50.0, 1.5]
for (name, est, true_val) in zip(config.theta_names, result.theta, true_values)
    bias = (est - true_val) / true_val * 100
    println("$(rpad(name, 12))$(lpad(string(round(est, digits=2)), 8))    $(true_val)   $(round(bias, digits=1))%")
end
