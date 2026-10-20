# Example 04: Coupled PK-PD Simulation
#
# This example demonstrates coupled pharmacokinetic-pharmacodynamic
# simulation using the indirect response turnover model.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "core", "NeoPKPD"))

using NeoPKPD

# === PK Model ===
# One-compartment IV bolus
pk_params = OneCompIVBolusParams(5.0, 50.0)  # CL=5, V=50
doses = [DoseEvent(0.0, 100.0)]
pk_spec = ModelSpec(OneCompIVBolus(), "pkpd_pk", pk_params, doses)

# === PD Model ===
# Indirect response with turnover
# Kin = 10 (input rate)
# Kout = 0.5 (output rate constant)
# R0 = 20 (baseline response = Kin/Kout)
# Imax = 0.8 (maximum inhibition)
# IC50 = 1.0 mg/L
pd_params = IndirectResponseTurnoverParams(10.0, 0.5, 20.0, 0.8, 1.0)
pd_spec = PDSpec(
    IndirectResponseTurnover(),
    "pkpd_pd",
    pd_params,
    :conc,      # PK observation to use as input
    :response   # PD output name
)

# Simulation over 72 hours (3 days) to see full PD response
grid = SimGrid(0.0, 72.0, collect(0.0:0.5:72.0))
solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10_000_000)

# Run coupled PKPD simulation
println("Running coupled PK-PD simulation...")
result = simulate_pkpd_coupled(pk_spec, pd_spec, grid, solver)

# Extract results
t = result.t
conc = result.observations[:conc]
response = result.observations[:response]

println("\n=== Coupled PK-PD Results ===")
println()
println("PK Parameters:")
println("  CL = 5.0 L/h, V = 50.0 L")
println("  t½ = $(round(0.693 * 50 / 5, digits=1)) h")
println()
println("PD Parameters:")
println("  Kin = 10, Kout = 0.5")
println("  R0 = 20 (baseline)")
println("  Imax = 0.8, IC50 = 1.0 mg/L")
println()

# Key metrics
c0 = conc[1]
r0 = response[1]
r_min = minimum(response)
r_min_idx = argmin(response)
t_min = t[r_min_idx]

println("Key Results:")
println("  Initial concentration: $(round(c0, digits=2)) mg/L")
println("  Initial response: $(round(r0, digits=2))")
println("  Minimum response: $(round(r_min, digits=2)) at t=$(round(t_min, digits=1)) h")
println("  Response suppression: $(round((r0 - r_min) / r0 * 100, digits=1))%")
println()

# Time course summary
println("Time Course:")
println("Time (h) | Conc (mg/L) | Response")
println("-" ^ 40)
for i in [1, 5, 9, 13, 25, 37, 49, 73, 97, 121, 145]
    if i <= length(t)
        println(
            "$(lpad(round(t[i], digits=0), 6)) | " *
            "$(lpad(round(conc[i], digits=3), 10)) | " *
            "$(round(response[i], digits=2))"
        )
    end
end

# Write artifact
outpath = joinpath(@__DIR__, "..", "output", "04_pkpd_coupled.json")
write_execution_json(
    outpath;
    model_spec=pk_spec,
    grid=grid,
    solver=solver,
    result=result,
    pd_spec=pd_spec
)
println("\nWrote: $outpath")
