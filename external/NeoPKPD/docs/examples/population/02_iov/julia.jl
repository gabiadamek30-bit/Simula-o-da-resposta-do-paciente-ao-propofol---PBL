# Inter-Occasion Variability (IOV) - Julia Example
# Run: julia --project=packages/core julia.jl

using NeoPKPD

println("Inter-Occasion Variability (IOV)")
println("="^50)

# Base model for single dose
base_model = create_model_spec(
    "OneCompOralFirstOrder";
    name = "iov_example",
    params = Dict("Ka" => 1.5, "CL" => 5.0, "V" => 50.0),
    doses = [DoseEvent(time=0.0, amount=100.0)]
)

# Population specification with IIV and IOV
n_subjects = 50
n_occasions = 3
pop_spec = create_population_spec(
    base_model_spec = base_model,
    iiv = LogNormalIIV(
        omegas = Dict("Ka" => 0.4, "CL" => 0.3, "V" => 0.2),
        seed = 12345
    ),
    iov = LogNormalIOV(
        pis = Dict("Ka" => 0.2, "CL" => 0.15),  # IOV on Ka and CL only
        occasions = n_occasions,
        seed = 54321
    ),
    n = n_subjects
)

# Simulation settings
grid = SimulationGrid(
    t0 = 0.0,
    t1 = 24.0,
    saveat = collect(0.0:0.5:24.0)
)

solver = SolverConfig(alg = "Tsit5", reltol = 1e-8, abstol = 1e-10)

# Run population simulation
println("\nRunning population simulation...")
println("  Subjects: $(n_subjects)")
println("  Occasions: $(n_occasions)")
result = simulate_population(pop_spec, grid, solver)

# Analyze variance components
println("\nVariance Components:")
println("-"^50)
println("Parameter   ω² (IIV)   π² (IOV)   Total Var")
println("-"^50)
println("Ka          0.16       0.04       0.20")
println("CL          0.09       0.0225     0.1125")
println("V           0.04       -          0.04")

println("\nConclusion: IOV adds to total variability within subjects")
println("Useful for modeling replicate PK studies")
