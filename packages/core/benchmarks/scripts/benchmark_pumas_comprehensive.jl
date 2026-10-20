#!/usr/bin/env julia
"""
Pumas Comprehensive Benchmark Suite
====================================

Prerequisites:
  - Pumas installed (requires commercial license from PumasAI)
  - Julia 1.9+

Installation:
  using Pkg
  Pkg.Registry.add(RegistrySpec(url="https://registry.pumas.ai"))
  Pkg.add("Pumas")
  Pkg.add("PumasUtilities")

Usage:
  julia benchmark_pumas_comprehensive.jl

Output:
  ../results/pumas_comprehensive_benchmarks.csv
"""

using Pkg

# Check for Pumas
try
    using Pumas
    using PumasUtilities
catch e
    println("ERROR: Pumas not found")
    println()
    println("Pumas is a commercial Julia package from PumasAI.")
    println("For installation (requires license):")
    println("  using Pkg")
    println("  Pkg.Registry.add(RegistrySpec(url=\"https://registry.pumas.ai\"))")
    println("  Pkg.add(\"Pumas\")")
    println()
    println("For benchmark comparison, use reference data from:")
    println("  - Pumas documentation: https://docs.pumas.ai/")
    println("  - Rackauckas C, et al. Pumas Performance. PAGE 2021")
    exit(1)
end

using Statistics
using Dates
using Printf
using Random

# =============================================================================
# Configuration
# =============================================================================

const SCRIPT_DIR = @__DIR__
const RESULTS_DIR = joinpath(SCRIPT_DIR, "..", "results")
const OUTPUT_FILE = joinpath(RESULTS_DIR, "pumas_comprehensive_benchmarks.csv")

const CONFIG = (
    n_runs = 50,
    n_warmup = 3,
    random_seed = 12345,
    timestamp = now(),
    pumas_version = string(pkgversion(Pumas)),
)

mkpath(RESULTS_DIR)

println("\n" * "="^70)
println("Pumas Comprehensive Benchmark Suite")
println("="^70)
println("Pumas Version: $(CONFIG.pumas_version)")
println("Julia Version: $(VERSION)")
println("Runs per benchmark: $(CONFIG.n_runs)")
println("Timestamp: $(CONFIG.timestamp)")
println("="^70)

# =============================================================================
# Benchmark Infrastructure
# =============================================================================

struct BenchmarkResult
    category::String
    subcategory::String
    name::String
    model::String
    n_subjects::Int
    times::Vector{Float64}
    timestamp::DateTime
    julia_version::String
    pumas_version::String
end

function summarize(r::BenchmarkResult)
    return (
        category = r.category,
        subcategory = r.subcategory,
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
        pumas_version = r.pumas_version,
    )
end

function run_benchmark(f::Function, name::String;
                       category::String, subcategory::String="",
                       model::String, n_subjects::Int=1)
    println("  Running: $name ($model)...")

    # Warmup
    for _ in 1:CONFIG.n_warmup
        try
            f()
        catch e
            println("    WARMUP ERROR: $e")
            return nothing
        end
    end

    # Timed runs
    times = Float64[]
    for _ in 1:CONFIG.n_runs
        GC.gc()
        t = @elapsed f()
        push!(times, t)
    end

    result = BenchmarkResult(
        category, subcategory, name, model, n_subjects, times,
        now(), string(VERSION), CONFIG.pumas_version
    )

    s = summarize(result)
    @printf("    Mean: %.3f ms (±%.3f), Median: %.3f ms\n",
            s.mean_ms, s.std_ms, s.median_ms)

    return result
end

# =============================================================================
# PK Models
# =============================================================================

function benchmark_pk_models()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PK MODELS")
    println("="^70)

    results = BenchmarkResult[]

    # One-Compartment IV Bolus
    model_1cmt = @model begin
        @param begin
            tvcl ∈ RealDomain(lower=0)
            tvv ∈ RealDomain(lower=0)
            Ω ∈ PDiagDomain(2)
            σ ∈ RealDomain(lower=0)
        end
        @random begin
            η ~ MvNormal(Ω)
        end
        @pre begin
            CL = tvcl * exp(η[1])
            V = tvv * exp(η[2])
        end
        @dynamics Central1
        @derived begin
            cp = @. Central / V
            dv ~ @. Normal(cp, σ)
        end
    end

    param_1cmt = (tvcl = 5.0, tvv = 50.0, Ω = Diagonal([0.09, 0.09]), σ = 0.1)

    dose = DosageRegimen(100, time = 0)
    subj = Subject(id = 1, events = dose)
    obstimes = 0:0.5:24

    r = run_benchmark(
        () -> simobs(model_1cmt, subj, param_1cmt, obstimes = obstimes),
        "simulate", category="pk", subcategory="compartmental",
        model="OneCompIVBolus"
    )
    r !== nothing && push!(results, r)

    # Two-Compartment IV Bolus
    model_2cmt = @model begin
        @param begin
            tvcl ∈ RealDomain(lower=0)
            tvv1 ∈ RealDomain(lower=0)
            tvq ∈ RealDomain(lower=0)
            tvv2 ∈ RealDomain(lower=0)
            Ω ∈ PDiagDomain(4)
            σ ∈ RealDomain(lower=0)
        end
        @random begin
            η ~ MvNormal(Ω)
        end
        @pre begin
            CL = tvcl * exp(η[1])
            V1 = tvv1 * exp(η[2])
            Q = tvq * exp(η[3])
            V2 = tvv2 * exp(η[4])
        end
        @dynamics Central1Periph1
        @derived begin
            cp = @. Central / V1
            dv ~ @. Normal(cp, σ)
        end
    end

    param_2cmt = (
        tvcl = 5.0, tvv1 = 50.0, tvq = 10.0, tvv2 = 100.0,
        Ω = Diagonal([0.09, 0.09, 0.09, 0.09]), σ = 0.1
    )

    r = run_benchmark(
        () -> simobs(model_2cmt, subj, param_2cmt, obstimes = obstimes),
        "simulate", category="pk", subcategory="compartmental",
        model="TwoCompIVBolus"
    )
    r !== nothing && push!(results, r)

    # One-Compartment Oral
    model_oral = @model begin
        @param begin
            tvka ∈ RealDomain(lower=0)
            tvcl ∈ RealDomain(lower=0)
            tvv ∈ RealDomain(lower=0)
            Ω ∈ PDiagDomain(3)
            σ ∈ RealDomain(lower=0)
        end
        @random begin
            η ~ MvNormal(Ω)
        end
        @pre begin
            Ka = tvka * exp(η[1])
            CL = tvcl * exp(η[2])
            V = tvv * exp(η[3])
        end
        @dynamics Depots1Central1
        @derived begin
            cp = @. Central / V
            dv ~ @. Normal(cp, σ)
        end
    end

    param_oral = (
        tvka = 1.5, tvcl = 5.0, tvv = 50.0,
        Ω = Diagonal([0.09, 0.09, 0.09]), σ = 0.1
    )

    r = run_benchmark(
        () -> simobs(model_oral, subj, param_oral, obstimes = obstimes),
        "simulate", category="pk", subcategory="compartmental",
        model="OneCompOralFirstOrder"
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# PD Models
# =============================================================================

function benchmark_pd_models()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PD MODELS")
    println("="^70)

    results = BenchmarkResult[]

    # Direct Emax model
    model_emax = @model begin
        @param begin
            e0 ∈ RealDomain()
            emax ∈ RealDomain(lower=0)
            ec50 ∈ RealDomain(lower=0)
            σ ∈ RealDomain(lower=0)
        end
        @pre begin
            E0 = e0
            Emax = emax
            EC50 = ec50
        end
        @vars begin
            conc = t  # Placeholder, would come from PK
        end
        @derived begin
            effect = @. E0 + Emax * conc / (EC50 + conc)
            dv ~ @. Normal(effect, σ)
        end
    end

    param_emax = (e0 = 0.0, emax = 100.0, ec50 = 10.0, σ = 0.1)

    subj = Subject(id = 1)
    obstimes = 0:0.5:24

    r = run_benchmark(
        () -> simobs(model_emax, subj, param_emax, obstimes = obstimes),
        "direct_effect", category="pd", subcategory="direct",
        model="DirectEmax"
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# PKPD Models
# =============================================================================

function benchmark_pkpd_models()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PKPD MODELS")
    println("="^70)

    results = BenchmarkResult[]

    # PK-PD linked model
    model_pkpd = @model begin
        @param begin
            tvcl ∈ RealDomain(lower=0)
            tvv ∈ RealDomain(lower=0)
            e0 ∈ RealDomain()
            emax ∈ RealDomain(lower=0)
            ec50 ∈ RealDomain(lower=0)
            Ω ∈ PDiagDomain(2)
            σ_pk ∈ RealDomain(lower=0)
            σ_pd ∈ RealDomain(lower=0)
        end
        @random begin
            η ~ MvNormal(Ω)
        end
        @pre begin
            CL = tvcl * exp(η[1])
            V = tvv * exp(η[2])
            E0 = e0
            Emax = emax
            EC50 = ec50
        end
        @dynamics Central1
        @derived begin
            cp = @. Central / V
            effect = @. E0 + Emax * cp / (EC50 + cp)
            dv_pk ~ @. Normal(cp, σ_pk)
            dv_pd ~ @. Normal(effect, σ_pd)
        end
    end

    param_pkpd = (
        tvcl = 5.0, tvv = 50.0,
        e0 = 0.0, emax = 100.0, ec50 = 10.0,
        Ω = Diagonal([0.09, 0.09]),
        σ_pk = 0.1, σ_pd = 0.1
    )

    dose = DosageRegimen(100, time = 0)
    subj = Subject(id = 1, events = dose)
    obstimes = 0:1.0:72

    r = run_benchmark(
        () -> simobs(model_pkpd, subj, param_pkpd, obstimes = obstimes),
        "pkpd_coupled", category="pkpd", subcategory="direct",
        model="OneComp_DirectEmax"
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# Population Simulation
# =============================================================================

function benchmark_population()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: POPULATION SIMULATION")
    println("="^70)

    results = BenchmarkResult[]

    model_pop = @model begin
        @param begin
            tvcl ∈ RealDomain(lower=0)
            tvv1 ∈ RealDomain(lower=0)
            tvq ∈ RealDomain(lower=0)
            tvv2 ∈ RealDomain(lower=0)
            Ω ∈ PDiagDomain(4)
            σ ∈ RealDomain(lower=0)
        end
        @random begin
            η ~ MvNormal(Ω)
        end
        @pre begin
            CL = tvcl * exp(η[1])
            V1 = tvv1 * exp(η[2])
            Q = tvq * exp(η[3])
            V2 = tvv2 * exp(η[4])
        end
        @dynamics Central1Periph1
        @derived begin
            cp = @. Central / V1
            dv ~ @. Normal(cp, σ)
        end
    end

    param_pop = (
        tvcl = 5.0, tvv1 = 50.0, tvq = 10.0, tvv2 = 100.0,
        Ω = Diagonal([0.09, 0.09, 0.09, 0.09]),
        σ = 0.1
    )

    obstimes = 0:0.5:24

    for n_subj in [10, 50, 100, 500, 1000]
        dose = DosageRegimen(100, time = 0)
        pop = [Subject(id = i, events = dose) for i in 1:n_subj]

        r = run_benchmark(
            () -> simobs(model_pop, pop, param_pop, obstimes = obstimes),
            "population_simulation", category="population", subcategory="iiv",
            model="TwoCompIV_IIV", n_subjects=n_subj
        )
        r !== nothing && push!(results, r)
    end

    return results
end

# =============================================================================
# NCA
# =============================================================================

function benchmark_nca()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: NCA")
    println("="^70)

    results = BenchmarkResult[]

    # Generate concentration-time data
    times = 0:0.5:24
    conc = 100 .* exp.(-0.1 .* times)

    # Create NCA dataset
    nca_data = NCASubject(conc, times, id = 1, dose = 100)

    r = run_benchmark(
        () -> NCA.auc(nca_data),
        "auc_0_t", category="nca", subcategory="exposure",
        model="AUC_0_t"
    )
    r !== nothing && push!(results, r)

    r = run_benchmark(
        () -> NCA.cmax(nca_data),
        "cmax_tmax", category="nca", subcategory="exposure",
        model="Cmax_Tmax"
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# Main Execution
# =============================================================================

function main()
    Random.seed!(CONFIG.random_seed)

    all_results = BenchmarkResult[]

    append!(all_results, benchmark_pk_models())
    append!(all_results, benchmark_pd_models())
    append!(all_results, benchmark_pkpd_models())
    append!(all_results, benchmark_population())
    append!(all_results, benchmark_nca())

    # Save results
    open(OUTPUT_FILE, "w") do io
        println(io, "\"category\",\"subcategory\",\"name\",\"model\",\"n_subjects\",\"n_runs\",\"mean_ms\",\"std_ms\",\"median_ms\",\"min_ms\",\"max_ms\",\"ci_lower_ms\",\"ci_upper_ms\",\"timestamp\",\"julia_version\",\"pumas_version\"")
        for r in all_results
            s = summarize(r)
            println(io, "\"$(s.category)\",\"$(s.subcategory)\",\"$(s.name)\",\"$(s.model)\",$(s.n_subjects),$(s.n_runs),$(s.mean_ms),$(s.std_ms),$(s.median_ms),$(s.min_ms),$(s.max_ms),$(s.ci_lower_ms),$(s.ci_upper_ms),\"$(s.timestamp)\",\"$(s.julia_version)\",\"$(s.pumas_version)\"")
        end
    end

    println("\n" * "="^70)
    println("PUMAS BENCHMARK COMPLETE")
    println("Results saved to: $OUTPUT_FILE")
    println("Total benchmarks: $(length(all_results))")
    println("="^70)
end

main()
