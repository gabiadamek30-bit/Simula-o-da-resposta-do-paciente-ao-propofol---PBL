#!/usr/bin/env julia
"""
NeoPKPD Benchmark Suite
=======================

Comprehensive benchmarks for NeoPKPD Julia core.
Produces reproducible timing data with statistical analysis.

Usage:
    julia --project=packages/core packages/core/benchmarks/scripts/benchmark_neopkpd.jl

Output:
    packages/core/benchmarks/results/neopkpd_julia_benchmarks.csv
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using NeoPKPD
using Statistics
using Dates
using Random
using Printf

# =============================================================================
# Configuration
# =============================================================================

const CONFIG = (
    n_runs = 100,           # Number of timed runs per benchmark
    n_warmup = 5,           # Warmup runs (JIT compilation)
    random_seed = 12345,    # For reproducibility
    output_dir = joinpath(@__DIR__, "..", "results"),
)

# =============================================================================
# Benchmark Infrastructure
# =============================================================================

struct BenchmarkResult
    category::String
    name::String
    model::String
    n_subjects::Int
    times::Vector{Float64}
    timestamp::DateTime
    julia_version::String
    neopkpd_version::String
end

function summarize(r::BenchmarkResult)
    return (
        category = r.category,
        name = r.name,
        model = r.model,
        n_subjects = r.n_subjects,
        n_runs = length(r.times),
        mean_ms = mean(r.times) * 1000,
        std_ms = std(r.times) * 1000,
        median_ms = median(r.times) * 1000,
        min_ms = minimum(r.times) * 1000,
        max_ms = maximum(r.times) * 1000,
        ci_lower_ms = quantile(r.times, 0.025) * 1000,
        ci_upper_ms = quantile(r.times, 0.975) * 1000,
        timestamp = r.timestamp,
        julia_version = r.julia_version,
        neopkpd_version = r.neopkpd_version,
    )
end

function run_benchmark(f::Function, name::String;
                       category::String, model::String, n_subjects::Int=1)
    println("  Running: $name ($model)...")

    # Warmup
    for _ in 1:CONFIG.n_warmup
        f()
    end

    # Timed runs
    times = Float64[]
    for i in 1:CONFIG.n_runs
        GC.gc()  # Clean GC before each run
        t = @elapsed f()
        push!(times, t)
    end

    result = BenchmarkResult(
        category, name, model, n_subjects, times,
        now(), string(VERSION), NEOPKPD_VERSION
    )

    s = summarize(result)
    @printf("    Mean: %.3f ms (±%.3f), Median: %.3f ms\n",
            s.mean_ms, s.std_ms, s.median_ms)

    return result
end

# =============================================================================
# Benchmark Definitions
# =============================================================================

function setup_common()
    return (
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0)),
        solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10_000_000),
        doses = [DoseEvent(0.0, 100.0)],
    )
end

# -----------------------------------------------------------------------------
# 1. Single Simulation Benchmarks
# -----------------------------------------------------------------------------

function benchmark_single_simulations()
    println("\n" * "="^60)
    println("BENCHMARK: Single Simulation Speed")
    println("="^60)

    common = setup_common()
    results = BenchmarkResult[]

    # 1a. One-compartment IV Bolus
    params_1comp = OneCompIVBolusParams(5.0, 50.0)  # CL=5, V=50
    spec_1comp = ModelSpec(OneCompIVBolus(), "bench_1comp", params_1comp, common.doses)

    push!(results, run_benchmark(
        () -> simulate(spec_1comp, common.grid, common.solver),
        "single_simulation", category="simulation", model="OneCompIV", n_subjects=1
    ))

    # 1b. Two-compartment IV Bolus
    params_2comp = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)  # CL, V1, Q, V2
    spec_2comp = ModelSpec(TwoCompIVBolus(), "bench_2comp", params_2comp, common.doses)

    push!(results, run_benchmark(
        () -> simulate(spec_2comp, common.grid, common.solver),
        "single_simulation", category="simulation", model="TwoCompIV", n_subjects=1
    ))

    # 1c. Two-compartment Oral
    params_2comp_oral = TwoCompOralParams(1.5, 5.0, 50.0, 10.0, 100.0)  # Ka, CL, V1, Q, V2
    spec_2comp_oral = ModelSpec(TwoCompOral(), "bench_2comp_oral", params_2comp_oral, common.doses)

    push!(results, run_benchmark(
        () -> simulate(spec_2comp_oral, common.grid, common.solver),
        "single_simulation", category="simulation", model="TwoCompOral", n_subjects=1
    ))

    # 1d. Transit Absorption
    params_transit = TransitAbsorptionParams(3, 0.5, 1.5, 5.0, 50.0)  # N, Ktr, Ka, CL, V
    spec_transit = ModelSpec(TransitAbsorption(), "bench_transit", params_transit, common.doses)

    push!(results, run_benchmark(
        () -> simulate(spec_transit, common.grid, common.solver),
        "single_simulation", category="simulation", model="TransitAbsorption", n_subjects=1
    ))

    # 1e. Michaelis-Menten
    params_mm = MichaelisMentenEliminationParams(50.0, 10.0, 50.0)  # Vmax, Km, V
    spec_mm = ModelSpec(MichaelisMentenElimination(), "bench_mm", params_mm, common.doses)

    push!(results, run_benchmark(
        () -> simulate(spec_mm, common.grid, common.solver),
        "single_simulation", category="simulation", model="MichaelisMenten", n_subjects=1
    ))

    return results
end

# -----------------------------------------------------------------------------
# 2. Population Simulation Benchmarks
# -----------------------------------------------------------------------------

function benchmark_population_simulations()
    println("\n" * "="^60)
    println("BENCHMARK: Population Simulation Speed")
    println("="^60)

    common = setup_common()
    results = BenchmarkResult[]

    # Base model for population
    params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)  # CL, V1, Q, V2
    spec = ModelSpec(TwoCompIVBolus(), "bench_pop", params, common.doses)

    for n_subj in [50, 100, 500, 1000]
        iiv = IIVSpec(
            LogNormalIIV(),
            Dict(:CL => 0.3, :V1 => 0.2),  # 30% CV for CL, 20% for V1
            UInt64(CONFIG.random_seed),
            n_subj
        )
        pop_spec = PopulationSpec(spec, iiv, nothing, nothing, IndividualCovariates[])

        push!(results, run_benchmark(
            () -> simulate_population(pop_spec, common.grid, common.solver),
            "population_simulation",
            category="population",
            model="TwoCompIV",
            n_subjects=n_subj
        ))
    end

    return results
end

# -----------------------------------------------------------------------------
# 3. PKPD Simulation Benchmarks
# -----------------------------------------------------------------------------

function benchmark_pkpd_simulations()
    println("\n" * "="^60)
    println("BENCHMARK: PKPD Coupled Simulation Speed")
    println("="^60)

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10_000_000)
    doses = [DoseEvent(0.0, 100.0)]
    results = BenchmarkResult[]

    # 3a. Direct Emax
    pk_params = OneCompIVBolusParams(5.0, 50.0)
    pd_params = DirectEmaxParams(0.0, 100.0, 10.0)  # E0, Emax, EC50

    push!(results, run_benchmark(
        () -> simulate_pkpd(
            OneCompIVBolus(), pk_params,
            DirectEmax(), pd_params,
            doses, grid, solver
        ),
        "pkpd_simulation", category="pkpd", model="DirectEmax", n_subjects=1
    ))

    # 3b. Sigmoid Emax
    pd_params_sig = SigmoidEmaxParams(0.0, 100.0, 10.0, 2.0)  # E0, Emax, EC50, gamma

    push!(results, run_benchmark(
        () -> simulate_pkpd(
            OneCompIVBolus(), pk_params,
            SigmoidEmax(), pd_params_sig,
            doses, grid, solver
        ),
        "pkpd_simulation", category="pkpd", model="SigmoidEmax", n_subjects=1
    ))

    return results
end

# -----------------------------------------------------------------------------
# 4. Sensitivity Analysis Benchmarks
# -----------------------------------------------------------------------------

function benchmark_sensitivity()
    println("\n" * "="^60)
    println("BENCHMARK: Sensitivity Analysis Speed")
    println("="^60)

    common = setup_common()
    results = BenchmarkResult[]

    params = OneCompIVBolusParams(5.0, 50.0)
    spec = ModelSpec(OneCompIVBolus(), "bench_sens", params, common.doses)

    # Local sensitivity (single parameter)
    pert = Perturbation(RelativePerturbation(), :CL, 0.1)
    plan = PerturbationPlan("local_sens", [pert])

    push!(results, run_benchmark(
        () -> run_sensitivity(spec, common.grid, common.solver; plan=plan, observation=:conc),
        "sensitivity_local", category="sensitivity", model="LocalSingle", n_subjects=1
    ))

    # Local sensitivity (multiple parameters)
    perts = [
        Perturbation(RelativePerturbation(), :CL, 0.1),
        Perturbation(RelativePerturbation(), :V, 0.1),
    ]
    plan_multi = PerturbationPlan("local_sens_multi", perts)

    push!(results, run_benchmark(
        () -> run_sensitivity(spec, common.grid, common.solver; plan=plan_multi, observation=:conc),
        "sensitivity_local_multi", category="sensitivity", model="LocalMulti", n_subjects=1
    ))

    return results
end

# =============================================================================
# Main Runner
# =============================================================================

function run_all_benchmarks()
    println("\n" * "="^70)
    println("NeoPKPD Benchmark Suite")
    println("="^70)
    println("Julia Version: ", VERSION)
    println("NeoPKPD Version: ", NEOPKPD_VERSION)
    println("Runs per benchmark: ", CONFIG.n_runs)
    println("Warmup runs: ", CONFIG.n_warmup)
    println("Random seed: ", CONFIG.random_seed)
    println("Timestamp: ", now())
    println("="^70)

    Random.seed!(CONFIG.random_seed)

    all_results = BenchmarkResult[]

    # Run all benchmark categories
    append!(all_results, benchmark_single_simulations())
    append!(all_results, benchmark_population_simulations())

    # Try PKPD benchmarks (may not be available in all versions)
    try
        append!(all_results, benchmark_pkpd_simulations())
    catch e
        println("\n  Skipping PKPD benchmarks: ", e)
    end

    # Try sensitivity benchmarks
    try
        append!(all_results, benchmark_sensitivity())
    catch e
        println("\n  Skipping sensitivity benchmarks: ", e)
    end

    # Save results
    mkpath(CONFIG.output_dir)
    output_file = joinpath(CONFIG.output_dir, "neopkpd_julia_benchmarks.csv")

    # Convert to table format
    summaries = [summarize(r) for r in all_results]

    open(output_file, "w") do io
        # Header
        println(io, "category,name,model,n_subjects,n_runs,mean_ms,std_ms,median_ms,min_ms,max_ms,ci_lower_ms,ci_upper_ms,timestamp,julia_version,neopkpd_version")

        # Data rows
        for s in summaries
            println(io, join([
                s.category, s.name, s.model, s.n_subjects, s.n_runs,
                @sprintf("%.4f", s.mean_ms),
                @sprintf("%.4f", s.std_ms),
                @sprintf("%.4f", s.median_ms),
                @sprintf("%.4f", s.min_ms),
                @sprintf("%.4f", s.max_ms),
                @sprintf("%.4f", s.ci_lower_ms),
                @sprintf("%.4f", s.ci_upper_ms),
                s.timestamp,
                s.julia_version,
                s.neopkpd_version
            ], ","))
        end
    end

    println("\n" * "="^70)
    println("BENCHMARK COMPLETE")
    println("Results saved to: ", output_file)
    println("="^70)

    # Print summary table
    println("\n" * "-"^70)
    println("SUMMARY")
    println("-"^70)
    @printf("%-25s %-20s %10s %10s\n", "Category", "Model", "Mean (ms)", "Std (ms)")
    println("-"^70)
    for s in summaries
        @printf("%-25s %-20s %10.3f %10.3f\n", s.category, s.model, s.mean_ms, s.std_ms)
    end
    println("-"^70)

    return summaries
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_all_benchmarks()
end
