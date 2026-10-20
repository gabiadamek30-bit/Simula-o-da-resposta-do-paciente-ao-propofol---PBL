# Single Subject Sensitivity - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Single Subject Sensitivity Analysis")
println("="^50)

# Base model
model = create_model_spec(
    "OneCompOralFirstOrder";
    params = Dict("Ka" => 1.5, "CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Compute sensitivity
result = compute_sensitivity(
    model;
    parameters = ["Ka", "CL", "V"],
    perturbation = 0.10,  # ±10%
    metrics = ["auc", "cmax", "tmax", "t_half"],
    t_end = 24.0
)

println("\nNominal Values:")
println("-"^40)
println("  Ka = $(model.params["Ka"]) 1/h")
println("  CL = $(model.params["CL"]) L/h")
println("  V  = $(model.params["V"]) L")

println("\nNominal Metrics:")
println("-"^40)
for (metric, value) in result.nominal_metrics
    println("  $(rpad(metric, 8)): $(round(value, digits=3))")
end

println("\nSensitivity Coefficients:")
println("-"^50)
println("Parameter    AUC      Cmax     Tmax     t_half")
println("-"^50)
for param in ["Ka", "CL", "V"]
    s = result.sensitivity[param]
    println("  $(rpad(param, 4))      $(lpad(string(round(s["auc"], digits=3)), 6))   " *
            "$(lpad(string(round(s["cmax"], digits=3)), 6))   " *
            "$(lpad(string(round(s["tmax"], digits=3)), 6))   " *
            "$(lpad(string(round(s["t_half"], digits=3)), 6))")
end

println("\nInterpretation:")
println("-"^50)
println("Coefficient = (% change in metric) / (% change in parameter)")
