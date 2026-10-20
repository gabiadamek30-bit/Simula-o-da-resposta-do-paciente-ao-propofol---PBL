# Basic NCA Analysis - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Basic Non-Compartmental Analysis (NCA)")
println("="^50)

# Sample concentration-time data (100 mg IV bolus)
times = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
concentrations = [2.0, 1.81, 1.64, 1.35, 0.91, 0.61, 0.41, 0.18, 0.03]
dose = 100.0  # mg

println("\nConcentration-Time Data:")
println("-"^30)
println("Time (h)   Conc (mg/L)")
println("-"^30)
for (t, c) in zip(times, concentrations)
    println("  $(lpad(string(t), 4))       $(round(c, digits=3))")
end

# Compute NCA metrics
result = compute_nca(
    times = times,
    concentrations = concentrations,
    dose = dose,
    route = :iv,  # IV administration
    method = :linear_log  # Linear up, log down trapezoidal
)

println("\n" * "="^50)
println("NCA RESULTS")
println("="^50)

println("\nExposure Metrics:")
println("-"^40)
println("Cmax      = $(round(result.cmax, digits=3)) mg/L")
println("Tmax      = $(round(result.tmax, digits=2)) h")
println("Clast     = $(round(result.clast, digits=4)) mg/L")
println("Tlast     = $(round(result.tlast, digits=1)) h")

println("\nAUC:")
println("-"^40)
println("AUC_0_t   = $(round(result.auc_0_t, digits=2)) mg·h/L")
println("AUC_0_inf = $(round(result.auc_0_inf, digits=2)) mg·h/L")
println("AUC_%ext  = $(round(result.auc_pct_extrapolated, digits=1))%")

println("\nElimination:")
println("-"^40)
println("λz        = $(round(result.lambda_z, digits=4)) 1/h")
println("t_half    = $(round(result.t_half, digits=2)) h")
println("R²        = $(round(result.r_squared, digits=4))")
println("N points  = $(result.n_points_lambda_z)")

println("\nClearance and Volume:")
println("-"^40)
println("CL        = $(round(result.cl, digits=2)) L/h")
println("Vz        = $(round(result.vz, digits=1)) L")
println("Vss       = $(round(result.vss, digits=1)) L")

println("\nMean Residence Time:")
println("-"^40)
println("MRT       = $(round(result.mrt, digits=2)) h")
println("AUMC_0_inf= $(round(result.aumc_0_inf, digits=1)) mg·h²/L")
