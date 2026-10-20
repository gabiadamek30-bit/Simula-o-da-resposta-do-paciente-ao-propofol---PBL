# Standard VPC - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Standard Visual Predictive Check (VPC)")
println("="^50)

# Load or simulate observed data (using theophylline-like data)
observed = create_observed_data(
    ids = repeat(1:12, inner=11),
    times = repeat([0, 0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24], 12),
    dv = [  # Simulated theophylline-like concentrations
        10.5, 9.2, 8.5, 7.3, 5.8, 4.2, 3.1, 2.2, 1.5, 0.9, 0.2,  # Subject 1
        8.2, 7.5, 7.0, 6.1, 4.9, 3.6, 2.7, 1.9, 1.3, 0.8, 0.15,  # Subject 2
        11.8, 10.5, 9.7, 8.3, 6.6, 4.8, 3.5, 2.5, 1.7, 1.0, 0.22, # etc.
        9.5, 8.6, 8.0, 6.9, 5.5, 4.0, 3.0, 2.1, 1.4, 0.85, 0.18,
        12.2, 10.8, 10.0, 8.6, 6.9, 5.0, 3.7, 2.6, 1.8, 1.1, 0.25,
        7.8, 7.1, 6.6, 5.7, 4.6, 3.3, 2.5, 1.8, 1.2, 0.72, 0.13,
        10.1, 9.1, 8.4, 7.2, 5.7, 4.1, 3.1, 2.2, 1.5, 0.88, 0.19,
        9.8, 8.9, 8.2, 7.1, 5.6, 4.1, 3.0, 2.1, 1.45, 0.87, 0.18,
        11.2, 10.0, 9.2, 7.9, 6.3, 4.6, 3.4, 2.4, 1.6, 0.98, 0.21,
        8.8, 8.0, 7.4, 6.4, 5.1, 3.7, 2.8, 1.95, 1.32, 0.79, 0.16,
        10.8, 9.7, 9.0, 7.7, 6.1, 4.4, 3.3, 2.3, 1.55, 0.94, 0.20,
        9.2, 8.3, 7.7, 6.6, 5.3, 3.8, 2.85, 2.0, 1.38, 0.82, 0.17
    ]
)

# Population model (estimated parameters)
base_model = create_model_spec(
    "OneCompOralFirstOrder";
    params = Dict("Ka" => 1.5, "CL" => 0.5, "V" => 30.0),
    doses = [DoseEvent(time=0.0, amount=320.0)]  # Theophylline dose
)

pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIV(
        omegas = Dict("Ka" => 0.4, "CL" => 0.25, "V" => 0.15),
        seed = 12345
    ),
    n = 12
)

# Generate VPC
println("\nGenerating VPC (500 simulations)...")
vpc_result = compute_vpc(
    observed,
    pop_spec,
    n_simulations = 500,
    quantiles = [0.05, 0.5, 0.95],
    seed = 12345
)

# Display results
println("\n" * "="^60)
println("VPC RESULTS")
println("="^60)

println("\nObserved Data Quantiles:")
println("-"^50)
println("Time (h)  5th%   Median  95th%")
println("-"^50)
for i in 1:length(vpc_result.time_bins)
    t = vpc_result.time_bins[i]
    obs_5 = vpc_result.observed_quantiles[1][i]
    obs_50 = vpc_result.observed_quantiles[2][i]
    obs_95 = vpc_result.observed_quantiles[3][i]
    println("  $(lpad(string(round(t, digits=1)), 4))   $(round(obs_5, digits=2))   $(round(obs_50, digits=2))    $(round(obs_95, digits=2))")
end

println("\nSimulated Prediction Intervals:")
println("-"^50)
println("Median CI: [$(round(minimum(vpc_result.simulated_ci_median[1]), digits=2)), $(round(maximum(vpc_result.simulated_ci_median[2]), digits=2))]")
println("5th% CI:   [$(round(minimum(vpc_result.simulated_ci_5th[1]), digits=2)), $(round(maximum(vpc_result.simulated_ci_5th[2]), digits=2))]")
println("95th% CI:  [$(round(minimum(vpc_result.simulated_ci_95th[1]), digits=2)), $(round(maximum(vpc_result.simulated_ci_95th[2]), digits=2))]")

println("\nVPC Diagnostics:")
println("-"^50)
pct_within_5 = sum(vpc_result.obs_within_sim_ci_5th) / length(vpc_result.obs_within_sim_ci_5th) * 100
pct_within_50 = sum(vpc_result.obs_within_sim_ci_median) / length(vpc_result.obs_within_sim_ci_median) * 100
pct_within_95 = sum(vpc_result.obs_within_sim_ci_95th) / length(vpc_result.obs_within_sim_ci_95th) * 100
println("% Obs within 5th% CI:     $(round(pct_within_5, digits=1))%")
println("% Obs within Median CI:   $(round(pct_within_50, digits=1))%")
println("% Obs within 95th% CI:    $(round(pct_within_95, digits=1))%")

println("\nConclusion:")
if pct_within_50 > 80 && pct_within_5 > 70 && pct_within_95 > 70
    println("  VPC indicates adequate model fit")
else
    println("  VPC suggests model misspecification - investigate further")
end
