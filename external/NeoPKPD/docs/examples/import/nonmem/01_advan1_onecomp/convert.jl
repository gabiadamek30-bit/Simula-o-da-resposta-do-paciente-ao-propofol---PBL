# NONMEM ADVAN1 Import - Julia Example
# Run: julia --project=packages/core convert.jl

using Pkg
Pkg.activate("packages/core")

using NeoPKPD
using JSON

println("NONMEM ADVAN1 Import")
println("="^50)

# Import NONMEM control file with dose
model = import_nonmem(
    joinpath(@__DIR__, "run001.ctl");
    doses=[DoseEvent(0.0, 100.0)]
)

println("\nImported Model:")
println("  Type:       $(model.model_type)")
println("  Parameters: $(model.parameters)")
println("  Source:     $(model.source)")

# Display NeoPKPD specification
println("\nNeoPKPD Specification:")
println("-"^50)
spec = model.spec

println("  Model:   $(typeof(spec.kind))")
println("  Params:  $(spec.params)")

if model.iiv !== nothing
    println("  IIV:")
    println("    Kind:   $(model.iiv.kind)")
    println("    Omegas: $(model.iiv.omegas)")
end

if model.error !== nothing
    println("  Error:   $(model.error)")
end

if !isempty(model.warnings)
    println("\n  Warnings:")
    for w in model.warnings
        println("    - $w")
    end
end

# Simulate to verify
println("\n" * "="^50)
println("SIMULATION TEST")
println("="^50)

grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)
result = simulate(spec, grid, solver)

println("\nSimulation successful: $(length(result.times)) time points")
println("Cmax: $(round(maximum(result.observations[:conc]), digits=3)) mg/L")

# Save results
output = Dict(
    "model_type" => String(model.model_type),
    "parameters" => Dict(String(k) => v for (k, v) in model.parameters),
    "source" => model.source,
    "cmax" => maximum(result.observations[:conc])
)
open(joinpath(@__DIR__, "output.json"), "w") do io
    JSON.print(io, output, 2)
end
println("\nWrote output.json")
