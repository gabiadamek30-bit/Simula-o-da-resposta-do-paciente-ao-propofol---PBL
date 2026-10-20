# Example 03: Population Simulation with IIV
#
# This example demonstrates population pharmacokinetic simulation
# with inter-individual variability (IIV).

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "core", "NeoPKPD"))

using NeoPKPD

# Base model parameters (typical values)
# CL = 5.0 L/h, V = 50.0 L
base_params = OneCompIVBolusParams(5.0, 50.0)

# 100 mg IV bolus
doses = [DoseEvent(0.0, 100.0)]

# Base model specification
base_spec = ModelSpec(OneCompIVBolus(), "pop_iiv_example", base_params, doses)

# Inter-Individual Variability (IIV) specification
# - 30% CV on CL (omega ≈ 0.3 for small CV)
# - 20% CV on V (omega ≈ 0.2)
# - Seed for reproducibility
# - 100 individuals
iiv = IIVSpec(
    LogNormalIIV(),
    Dict(:CL => 0.3, :V => 0.2),
    UInt64(12345),
    100
)

# Create population specification
pop = PopulationSpec(base_spec, iiv, nothing, nothing, [])

# Simulation grid
grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run population simulation
println("Running population simulation with 100 individuals...")
result = simulate_population(pop, grid, solver)

# Analyze results
summary = result.summaries[:conc]

println("\n=== Population PK Results ===")
println("Number of individuals: $(length(result.individuals))")
println()

# Parameter distribution
cls = [p[:CL] for p in result.params]
vs = [p[:V] for p in result.params]
println("CL distribution:")
println("  Mean: $(round(sum(cls)/length(cls), digits=2)) L/h")
println("  Range: $(round(minimum(cls), digits=2)) - $(round(maximum(cls), digits=2)) L/h")
println()
println("V distribution:")
println("  Mean: $(round(sum(vs)/length(vs), digits=2)) L")
println("  Range: $(round(minimum(vs), digits=2)) - $(round(maximum(vs), digits=2)) L")
println()

# Concentration summary at selected time points
println("Concentration summary (mg/L):")
println("Time (h) | Mean    | Median  | 5%      | 95%")
println("-" ^ 50)
for i in [1, 5, 9, 13, 17, 21, 25, 33, 41, 49]
    if i <= length(summary.mean)
        t_val = result.individuals[1].t[i]
        println(
            "$(lpad(round(t_val, digits=1), 6))   | " *
            "$(lpad(round(summary.mean[i], digits=3), 7)) | " *
            "$(lpad(round(summary.median[i], digits=3), 7)) | " *
            "$(lpad(round(summary.quantiles[0.05][i], digits=3), 7)) | " *
            "$(lpad(round(summary.quantiles[0.95][i], digits=3), 7))"
        )
    end
end

# Write artifact
outpath = joinpath(@__DIR__, "..", "output", "03_population_iiv.json")
write_population_json(
    outpath;
    population_spec=pop,
    grid=grid,
    solver=solver,
    result=result
)
println("\nWrote: $outpath")
