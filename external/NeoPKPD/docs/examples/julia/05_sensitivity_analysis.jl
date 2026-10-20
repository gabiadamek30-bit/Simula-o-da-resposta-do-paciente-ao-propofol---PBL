# Example 05: Parameter Sensitivity Analysis
#
# This example demonstrates sensitivity analysis to understand
# how parameter changes affect model outputs.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "core", "NeoPKPD"))

using NeoPKPD

# Base model
params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)  # Ka, CL, V
doses = [DoseEvent(0.0, 200.0)]
spec = ModelSpec(OneCompOralFirstOrder(), "sensitivity_example", params, doses)

grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Define perturbation plans for each parameter
plans = [
    PerturbationPlan("Ka+10%", [Perturbation(RelativePerturbation(), :Ka, 0.1)]),
    PerturbationPlan("Ka-10%", [Perturbation(RelativePerturbation(), :Ka, -0.1)]),
    PerturbationPlan("CL+10%", [Perturbation(RelativePerturbation(), :CL, 0.1)]),
    PerturbationPlan("CL-10%", [Perturbation(RelativePerturbation(), :CL, -0.1)]),
    PerturbationPlan("V+10%", [Perturbation(RelativePerturbation(), :V, 0.1)]),
    PerturbationPlan("V-10%", [Perturbation(RelativePerturbation(), :V, -0.1)]),
]

println("=== Parameter Sensitivity Analysis ===")
println()
println("Base parameters: Ka=1.5/h, CL=5.0 L/h, V=50.0 L")
println("Perturbation: ±10%")
println()

# Run sensitivity analysis for each plan
results = Dict{String, SensitivityResult}()
for plan in plans
    results[plan.name] = run_sensitivity(spec, grid, solver; plan=plan, observation=:conc)
end

# Display results
println("Sensitivity Results (on concentration):")
println("=" ^ 60)
println("Plan        | Max Abs Δ | Max Rel Δ | L2 Norm Δ")
println("-" ^ 60)

for plan in plans
    r = results[plan.name]
    m = r.metrics
    println(
        "$(rpad(plan.name, 11)) | " *
        "$(lpad(round(m.max_abs_delta, digits=4), 9)) | " *
        "$(lpad(round(m.max_rel_delta, digits=4), 9)) | " *
        "$(lpad(round(m.l2_norm_delta, digits=4), 9))"
    )
end

# Rank parameters by sensitivity
println()
println("Parameter Sensitivity Ranking (by max relative delta):")
println("-" ^ 40)

# Get max sensitivity for each parameter (average of +/- perturbations)
param_sens = Dict{Symbol, Float64}()
for sym in [:Ka, :CL, :V]
    plus_plan = "$(sym)+10%"
    minus_plan = "$(sym)-10%"
    avg_sens = (results[plus_plan].metrics.max_rel_delta +
                results[minus_plan].metrics.max_rel_delta) / 2
    param_sens[sym] = avg_sens
end

sorted_params = sort(collect(param_sens), by=x->x[2], rev=true)
for (i, (param, sens)) in enumerate(sorted_params)
    println("$i. $param: $(round(sens * 100, digits=2))% max relative change")
end

# Write one sensitivity artifact as example
outpath = joinpath(@__DIR__, "..", "output", "05_sensitivity_CL.json")
write_sensitivity_json(
    outpath;
    model_spec=spec,
    grid=grid,
    solver=solver,
    result=results["CL+10%"]
)
println("\nWrote: $outpath")
