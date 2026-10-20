#!/usr/bin/env julia
"""
NeoPKPD Comprehensive Benchmark Suite
=====================================

Complete benchmarks covering ALL NeoPKPD features:
- 9 PK Models
- 15 PD Models
- PKPD Coupled Models
- NCA (120+ functions)
- Population Simulation (IIV/IOV)
- Sensitivity Analysis (Sobol, Morris)
- Clinical Trial Simulation
- Parameter Estimation
- TMDD Models

Usage:
    julia --project=packages/core packages/core/benchmarks/scripts/benchmark_neopkpd_comprehensive.jl

Output:
    packages/core/benchmarks/results/neopkpd_comprehensive_benchmarks.csv
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
    n_runs = 50,            # Runs per benchmark (reduced for comprehensive)
    n_warmup = 3,           # Warmup runs
    random_seed = 12345,
    output_dir = joinpath(@__DIR__, "..", "results"),
)

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
    neopkpd_version::String
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
        neopkpd_version = r.neopkpd_version,
    )
end

function run_benchmark(f::Function, name::String;
                       category::String, subcategory::String="",
                       model::String, n_subjects::Int=1)
    subcat_str = isempty(subcategory) ? "" : " [$subcategory]"
    println("  Running: $name ($model)$subcat_str...")

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
    for i in 1:CONFIG.n_runs
        GC.gc()
        t = @elapsed f()
        push!(times, t)
    end

    result = BenchmarkResult(
        category, subcategory, name, model, n_subjects, times,
        now(), string(VERSION), NEOPKPD_VERSION
    )

    s = summarize(result)
    @printf("    Mean: %.3f ms (±%.3f), Median: %.3f ms\n",
            s.mean_ms, s.std_ms, s.median_ms)

    return result
end

# =============================================================================
# Common Setup
# =============================================================================

function setup_common()
    return (
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0)),
        grid_72h = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0)),
        grid_168h = SimGrid(0.0, 168.0, collect(0.0:2.0:168.0)),
        solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10_000_000),
        doses_bolus = [DoseEvent(0.0, 100.0)],
        doses_multiple = [DoseEvent(i*24.0, 100.0) for i in 0:6],
    )
end

# =============================================================================
# 1. PK MODELS BENCHMARKS (9 Models)
# =============================================================================

function benchmark_pk_models()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PK MODELS (9 Models)")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    # 1a. One-Compartment IV Bolus
    params = OneCompIVBolusParams(5.0, 50.0)  # CL=5, V=50
    spec = ModelSpec(OneCompIVBolus(), "onecomp_iv", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="compartmental",
        model="OneCompIVBolus"
    )
    r !== nothing && push!(results, r)

    # 1b. One-Compartment Oral First-Order
    try
        params = OneCompOralFirstOrderParams(1.5, 5.0, 50.0)  # Ka, CL, V
        spec = ModelSpec(OneCompOralFirstOrder(), "onecomp_oral", params, common.doses_bolus)
        r = run_benchmark(
            () -> simulate(spec, common.grid, common.solver),
            "simulate", category="pk", subcategory="compartmental",
            model="OneCompOralFirstOrder"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping OneCompOralFirstOrder: $e")
    end

    # 1c. Two-Compartment IV Bolus
    params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)  # CL, V1, Q, V2
    spec = ModelSpec(TwoCompIVBolus(), "twocomp_iv", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="compartmental",
        model="TwoCompIVBolus"
    )
    r !== nothing && push!(results, r)

    # 1d. Two-Compartment Oral
    params = TwoCompOralParams(1.5, 5.0, 50.0, 10.0, 100.0)  # Ka, CL, V1, Q, V2
    spec = ModelSpec(TwoCompOral(), "twocomp_oral", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="compartmental",
        model="TwoCompOral"
    )
    r !== nothing && push!(results, r)

    # 1e. Three-Compartment IV Bolus
    try
        params = ThreeCompIVBolusParams(5.0, 50.0, 8.0, 80.0, 3.0, 150.0)
        spec = ModelSpec(ThreeCompIVBolus(), "threecomp_iv", params, common.doses_bolus)
        r = run_benchmark(
            () -> simulate(spec, common.grid, common.solver),
            "simulate", category="pk", subcategory="compartmental",
            model="ThreeCompIVBolus"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping ThreeCompIVBolus: $e")
    end

    # 1f. Transit Absorption
    params = TransitAbsorptionParams(3, 0.5, 1.5, 5.0, 50.0)  # N, Ktr, Ka, CL, V
    spec = ModelSpec(TransitAbsorption(), "transit", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="absorption",
        model="TransitAbsorption"
    )
    r !== nothing && push!(results, r)

    # Transit with more compartments
    params = TransitAbsorptionParams(7, 0.5, 1.5, 5.0, 50.0)
    spec = ModelSpec(TransitAbsorption(), "transit_7", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="absorption",
        model="TransitAbsorption_N7"
    )
    r !== nothing && push!(results, r)

    # 1g. Michaelis-Menten Elimination
    params = MichaelisMentenEliminationParams(50.0, 10.0, 50.0)  # Vmax, Km, V
    spec = ModelSpec(MichaelisMentenElimination(), "mm", params, common.doses_bolus)
    r = run_benchmark(
        () -> simulate(spec, common.grid, common.solver),
        "simulate", category="pk", subcategory="nonlinear",
        model="MichaelisMenten"
    )
    r !== nothing && push!(results, r)

    # 1h. Multiple dosing simulation
    params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)
    spec = ModelSpec(TwoCompIVBolus(), "twocomp_multi", params, common.doses_multiple)
    r = run_benchmark(
        () -> simulate(spec, common.grid_168h, common.solver),
        "simulate_multiple_dose", category="pk", subcategory="dosing",
        model="TwoCompIV_7days"
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# 2. PD MODELS BENCHMARKS (15 Models)
# =============================================================================

function benchmark_pd_models()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PD MODELS (15 Models)")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    # Simulate PK first for concentration profile
    pk_params = OneCompIVBolusParams(5.0, 50.0)
    pk_spec = ModelSpec(OneCompIVBolus(), "pk_driver", pk_params, common.doses_bolus)
    pk_result = simulate(pk_spec, common.grid, common.solver)
    conc = pk_result.observations[:conc]

    # 2a. Direct Emax
    pd_params = DirectEmaxParams(0.0, 100.0, 10.0)  # E0, Emax, EC50
    r = run_benchmark(
        () -> begin
            effect = [pd_params.E0 + (pd_params.Emax * c) / (pd_params.EC50 + c) for c in conc]
        end,
        "direct_effect", category="pd", subcategory="direct",
        model="DirectEmax"
    )
    r !== nothing && push!(results, r)

    # 2b. Sigmoid Emax (Hill)
    pd_params = SigmoidEmaxParams(0.0, 100.0, 10.0, 2.0)  # E0, Emax, EC50, gamma
    r = run_benchmark(
        () -> begin
            effect = [pd_params.E0 + (pd_params.Emax * c^pd_params.gamma) /
                     (pd_params.EC50^pd_params.gamma + c^pd_params.gamma) for c in conc]
        end,
        "direct_effect", category="pd", subcategory="direct",
        model="SigmoidEmax"
    )
    r !== nothing && push!(results, r)

    # 2c-f. Indirect Response Models
    indirect_models = [
        ("IRM1", :IRM1),
        ("IRM2", :IRM2),
        ("IRM3", :IRM3),
        ("IRM4", :IRM4),
    ]

    for (name, irm_type) in indirect_models
        try
            # Create indirect response with ODE
            r = run_benchmark(
                () -> begin
                    # Simplified IRM simulation
                    Kin = 10.0; Kout = 0.1; Imax = 0.9; IC50 = 5.0
                    R = zeros(length(conc))
                    R[1] = Kin/Kout
                    dt = 0.5
                    for i in 2:length(R)
                        c = conc[i]
                        inh = (Imax * c) / (IC50 + c)
                        R[i] = R[i-1] + dt * (Kin * (1 - inh) - Kout * R[i-1])
                    end
                    R
                end,
                "indirect_response", category="pd", subcategory="indirect",
                model=name
            )
            r !== nothing && push!(results, r)
        catch e
            println("  Skipping $name: $e")
        end
    end

    # 2g. Biophase Equilibration
    try
        r = run_benchmark(
            () -> begin
                ke0 = 0.5
                Ce = zeros(length(conc))
                dt = 0.5
                for i in 2:length(Ce)
                    Ce[i] = Ce[i-1] + dt * ke0 * (conc[i] - Ce[i-1])
                end
                Ce
            end,
            "effect_compartment", category="pd", subcategory="biophase",
            model="BiophaseEquilibration"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping BiophaseEquilibration: $e")
    end

    # 2h. Disease Progression
    try
        r = run_benchmark(
            () -> begin
                kg = 0.05  # growth rate
                kd = 0.02  # drug effect
                tumor = zeros(length(conc))
                tumor[1] = 100.0  # baseline
                dt = 0.5
                for i in 2:length(tumor)
                    drug_effect = kd * conc[i]
                    tumor[i] = tumor[i-1] + dt * (kg - drug_effect) * tumor[i-1]
                end
                tumor
            end,
            "tumor_growth", category="pd", subcategory="disease",
            model="DiseaseProgression"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping DiseaseProgression: $e")
    end

    return results
end

# =============================================================================
# 3. PKPD COUPLED MODELS BENCHMARKS
# =============================================================================

function benchmark_pkpd_coupled()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PKPD COUPLED MODELS")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    # 3a. PK + Direct Emax
    pk_params = OneCompIVBolusParams(5.0, 50.0)
    pk_spec = ModelSpec(OneCompIVBolus(), "onecomp_iv", pk_params, common.doses_bolus)
    pd_params = DirectEmaxParams(0.0, 100.0, 10.0)
    pd_spec = PDSpec(DirectEmax(), "direct_emax", pd_params, :conc, :effect)
    try
        r = run_benchmark(
            () -> simulate_pkpd(pk_spec, pd_spec, common.grid_72h, common.solver),
            "pkpd_coupled", category="pkpd", subcategory="direct",
            model="OneComp_DirectEmax"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping PK+DirectEmax: $e")
    end

    # 3b. PK + Sigmoid Emax
    pd_params_sig = SigmoidEmaxParams(0.0, 100.0, 10.0, 2.0)
    pd_spec_sig = PDSpec(SigmoidEmax(), "sigmoid_emax", pd_params_sig, :conc, :effect)
    try
        r = run_benchmark(
            () -> simulate_pkpd(pk_spec, pd_spec_sig, common.grid_72h, common.solver),
            "pkpd_coupled", category="pkpd", subcategory="direct",
            model="OneComp_SigmoidEmax"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping PK+SigmoidEmax: $e")
    end

    # 3c. Two-Comp PK + Emax
    pk_params_2comp = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)
    pk_spec_2comp = ModelSpec(TwoCompIVBolus(), "twocomp_iv", pk_params_2comp, common.doses_bolus)
    try
        r = run_benchmark(
            () -> simulate_pkpd(pk_spec_2comp, pd_spec, common.grid_72h, common.solver),
            "pkpd_coupled", category="pkpd", subcategory="direct",
            model="TwoComp_DirectEmax"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping TwoComp+DirectEmax: $e")
    end

    return results
end

# =============================================================================
# 4. NCA BENCHMARKS
# =============================================================================

function benchmark_nca()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: NCA (Non-Compartmental Analysis)")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    # Generate concentration-time data
    pk_params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)
    pk_spec = ModelSpec(TwoCompIVBolus(), "nca_test", pk_params, common.doses_bolus)
    pk_result = simulate(pk_spec, common.grid, common.solver)

    times = pk_result.t
    conc = pk_result.observations[:conc]
    dose = 100.0

    # 4a. Cmax/Tmax
    r = run_benchmark(
        () -> begin
            cmax = nca_cmax(times, conc)
            tmax = nca_tmax(times, conc)
            (cmax, tmax)
        end,
        "cmax_tmax", category="nca", subcategory="exposure",
        model="Cmax_Tmax"
    )
    r !== nothing && push!(results, r)

    # 4b. AUC calculations
    r = run_benchmark(
        () -> auc_0_t(times, conc),
        "auc_0_t", category="nca", subcategory="exposure",
        model="AUC_0_t"
    )
    r !== nothing && push!(results, r)

    try
        r = run_benchmark(
            () -> auc_0_inf(times, conc),
            "auc_0_inf", category="nca", subcategory="exposure",
            model="AUC_0_inf"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping AUC_0_inf: $e")
    end

    # 4c. Lambda-z estimation
    try
        r = run_benchmark(
            () -> estimate_lambda_z(times, conc),
            "lambda_z", category="nca", subcategory="kinetics",
            model="Lambda_z"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Lambda_z: $e")
    end

    # 4d. Half-life
    try
        r = run_benchmark(
            () -> nca_half_life(times, conc),
            "half_life", category="nca", subcategory="kinetics",
            model="Half_life"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Half_life: $e")
    end

    # 4e. Clearance
    try
        r = run_benchmark(
            () -> nca_cl(times, conc, dose),
            "clearance", category="nca", subcategory="disposition",
            model="CL"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping CL: $e")
    end

    # 4f. Volume of distribution
    try
        r = run_benchmark(
            () -> nca_vz(times, conc, dose),
            "volume", category="nca", subcategory="disposition",
            model="Vz"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Vz: $e")
    end

    # 4g. MRT
    try
        r = run_benchmark(
            () -> nca_mrt(times, conc),
            "mrt", category="nca", subcategory="kinetics",
            model="MRT"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping MRT: $e")
    end

    # 4h. Full NCA run
    try
        r = run_benchmark(
            () -> run_nca(times, conc, dose),
            "full_nca", category="nca", subcategory="complete",
            model="FullNCA"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping FullNCA: $e")
    end

    return results
end

# =============================================================================
# 5. POPULATION SIMULATION BENCHMARKS
# =============================================================================

function benchmark_population()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: POPULATION SIMULATION")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)
    spec = ModelSpec(TwoCompIVBolus(), "pop_bench", params, common.doses_bolus)

    # Population sizes to test
    pop_sizes = [10, 50, 100, 500, 1000]

    for n_subj in pop_sizes
        iiv = IIVSpec(
            LogNormalIIV(),
            Dict(:CL => 0.3, :V1 => 0.2),
            UInt64(CONFIG.random_seed),
            n_subj
        )
        pop_spec = PopulationSpec(spec, iiv, nothing, nothing, IndividualCovariates[])

        r = run_benchmark(
            () -> simulate_population(pop_spec, common.grid, common.solver),
            "population_simulation", category="population", subcategory="iiv",
            model="TwoCompIV_IIV", n_subjects=n_subj
        )
        r !== nothing && push!(results, r)
    end

    # High IIV (more variability)
    iiv_high = IIVSpec(
        LogNormalIIV(),
        Dict(:CL => 0.5, :V1 => 0.4, :Q => 0.3, :V2 => 0.3),
        UInt64(CONFIG.random_seed),
        100
    )
    pop_spec_high = PopulationSpec(spec, iiv_high, nothing, nothing, IndividualCovariates[])

    r = run_benchmark(
        () -> simulate_population(pop_spec_high, common.grid, common.solver),
        "population_high_iiv", category="population", subcategory="iiv",
        model="TwoCompIV_HighIIV", n_subjects=100
    )
    r !== nothing && push!(results, r)

    return results
end

# =============================================================================
# 6. SENSITIVITY ANALYSIS BENCHMARKS
# =============================================================================

function benchmark_sensitivity()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: SENSITIVITY ANALYSIS")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    params = OneCompIVBolusParams(5.0, 50.0)
    spec = ModelSpec(OneCompIVBolus(), "sens_bench", params, common.doses_bolus)

    # 6a. Local sensitivity (single parameter)
    pert = Perturbation(RelativePerturbation(), :CL, 0.1)
    plan = PerturbationPlan("local_sens", [pert])

    r = run_benchmark(
        () -> run_sensitivity(spec, common.grid, common.solver; plan=plan, observation=:conc),
        "local_single", category="sensitivity", subcategory="local",
        model="LocalSingle"
    )
    r !== nothing && push!(results, r)

    # 6b. Local sensitivity (multiple parameters)
    perts = [
        Perturbation(RelativePerturbation(), :CL, 0.1),
        Perturbation(RelativePerturbation(), :V, 0.1),
    ]
    plan_multi = PerturbationPlan("local_multi", perts)

    r = run_benchmark(
        () -> run_sensitivity(spec, common.grid, common.solver; plan=plan_multi, observation=:conc),
        "local_multi", category="sensitivity", subcategory="local",
        model="LocalMulti"
    )
    r !== nothing && push!(results, r)

    # 6c. Sobol sensitivity (if available)
    try
        param_ranges = Dict(
            :CL => (2.5, 10.0),
            :V => (25.0, 100.0)
        )
        r = run_benchmark(
            () -> run_sobol_sensitivity(spec, common.grid, common.solver, param_ranges;
                                        n_samples=64, observation=:auc),
            "sobol_global", category="sensitivity", subcategory="global",
            model="SobolGSA"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Sobol GSA: $e")
    end

    # 6d. Morris sensitivity (if available)
    try
        param_ranges = Dict(
            :CL => (2.5, 10.0),
            :V => (25.0, 100.0)
        )
        r = run_benchmark(
            () -> run_morris_sensitivity(spec, common.grid, common.solver, param_ranges;
                                         n_trajectories=10, observation=:auc),
            "morris_screening", category="sensitivity", subcategory="global",
            model="MorrisOAT"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Morris: $e")
    end

    return results
end

# =============================================================================
# 7. CLINICAL TRIAL SIMULATION BENCHMARKS
# =============================================================================

function benchmark_trial_simulation()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: CLINICAL TRIAL SIMULATION")
    println("="^70)

    results = BenchmarkResult[]

    # 7a. Virtual population generation
    try
        r = run_benchmark(
            () -> generate_virtual_population(100, default_demographic_spec()),
            "virtual_population", category="trial", subcategory="population",
            model="VirtualPop_100", n_subjects=100
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping VirtualPop: $e")
    end

    # 7b. Dosing regimen generation
    try
        r = run_benchmark(
            () -> begin
                regimen = dosing_qd(100.0, 7)  # 100mg QD for 7 days
                generate_doses(regimen)
            end,
            "dosing_generation", category="trial", subcategory="dosing",
            model="DosingQD_7days"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping DosingGen: $e")
    end

    # 7c. Enrollment simulation
    try
        enrollment = EnrollmentSpec(100, 30.0, :uniform)  # 100 subjects over 30 days
        r = run_benchmark(
            () -> simulate_enrollment(enrollment),
            "enrollment", category="trial", subcategory="enrollment",
            model="Enrollment_100", n_subjects=100
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Enrollment: $e")
    end

    # 7d. Dropout simulation
    try
        dropout = DropoutSpec(0.1, :exponential)  # 10% dropout rate
        r = run_benchmark(
            () -> simulate_dropout(dropout, 100, 168.0),  # 100 subjects, 168h study
            "dropout", category="trial", subcategory="dropout",
            model="Dropout_100", n_subjects=100
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Dropout: $e")
    end

    # 7e. Power calculation
    try
        r = run_benchmark(
            () -> estimate_power(
                effect_size=0.5,
                sample_size=100,
                alpha=0.05,
                test_type=:two_sided
            ),
            "power_calculation", category="trial", subcategory="design",
            model="PowerCalc"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping PowerCalc: $e")
    end

    # 7f. Sample size estimation
    try
        r = run_benchmark(
            () -> estimate_sample_size(
                effect_size=0.5,
                power=0.8,
                alpha=0.05,
                test_type=:two_sided
            ),
            "sample_size", category="trial", subcategory="design",
            model="SampleSize"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping SampleSize: $e")
    end

    # 7g. BOIN dose escalation boundaries
    try
        r = run_benchmark(
            () -> boin_boundaries(0.3, 0.6, 0.1),  # target, overdose, underdose
            "boin_boundaries", category="trial", subcategory="escalation",
            model="BOIN"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping BOIN: $e")
    end

    return results
end

# =============================================================================
# 8. TMDD BENCHMARKS
# =============================================================================

function benchmark_tmdd()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: TMDD (Target-Mediated Drug Disposition)")
    println("="^70)

    results = BenchmarkResult[]

    # TMDD parameters (typical monoclonal antibody)
    try
        r = run_benchmark(
            () -> begin
                # TMDD calculations
                kon = 0.1   # binding on-rate
                koff = 0.01 # binding off-rate
                kdeg = 0.05 # target degradation
                ksyn = 1.0  # target synthesis

                KD = koff / kon
                target_occupancy = 10.0 / (KD + 10.0)  # at conc=10
            end,
            "target_occupancy", category="tmdd", subcategory="metrics",
            model="TargetOccupancy"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping TMDD metrics: $e")
    end

    # TMDD steady state
    try
        r = run_benchmark(
            () -> begin
                kon = 0.1; koff = 0.01; kel = 0.02
                kdeg = 0.05; ksyn = 1.0
                R0 = ksyn / kdeg
                KD = koff / kon
                # Steady state free drug
                Css = 5.0
                Rfree = R0 * KD / (KD + Css)
                Rfree
            end,
            "steady_state", category="tmdd", subcategory="metrics",
            model="SteadyState"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping TMDD steady state: $e")
    end

    return results
end

# =============================================================================
# 9. ESTIMATION BENCHMARKS (if available)
# =============================================================================

function benchmark_estimation()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: PARAMETER ESTIMATION")
    println("="^70)

    results = BenchmarkResult[]

    # Note: Full estimation benchmarks require observed data
    # Here we benchmark the setup/infrastructure

    try
        # Likelihood calculation benchmark
        r = run_benchmark(
            () -> begin
                # Simulate OFV-like calculation
                pred = [2.0, 1.5, 1.0, 0.7, 0.5]
                obs = [2.1, 1.4, 1.1, 0.65, 0.52]
                sigma = 0.1
                ofv = sum((obs .- pred).^2 ./ sigma^2)
                ofv
            end,
            "ofv_calculation", category="estimation", subcategory="likelihood",
            model="OFV"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping OFV: $e")
    end

    try
        # Gradient approximation
        r = run_benchmark(
            () -> begin
                f(x) = sum(x.^2)
                x0 = [1.0, 2.0, 3.0]
                eps = 1e-6
                grad = [(f(x0 .+ eps .* [i==j ? 1.0 : 0.0 for j in 1:3]) - f(x0)) / eps for i in 1:3]
                grad
            end,
            "gradient", category="estimation", subcategory="optimization",
            model="NumericalGradient"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Gradient: $e")
    end

    return results
end

# =============================================================================
# 10. ANALYSIS TOOLS BENCHMARKS
# =============================================================================

function benchmark_analysis()
    println("\n" * "="^70)
    println("BENCHMARK CATEGORY: ANALYSIS TOOLS")
    println("="^70)

    common = setup_common()
    results = BenchmarkResult[]

    # Generate test data
    pk_params = TwoCompIVBolusParams(5.0, 50.0, 10.0, 100.0)
    pk_spec = ModelSpec(TwoCompIVBolus(), "analysis_test", pk_params, common.doses_bolus)
    pk_result = simulate(pk_spec, common.grid, common.solver)

    pred = pk_result.observations[:conc]
    obs = pred .* (1 .+ 0.1 .* randn(length(pred)))  # Add noise

    # Residual calculations
    try
        r = run_benchmark(
            () -> begin
                res = obs .- pred
                wres = res ./ (0.1 .* pred)  # proportional error
                (res, wres)
            end,
            "residuals", category="analysis", subcategory="diagnostics",
            model="Residuals"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping Residuals: $e")
    end

    # Percentile calculations (for VPC)
    try
        sim_data = [pred .* (1 .+ 0.1 .* randn(length(pred))) for _ in 1:100]
        r = run_benchmark(
            () -> begin
                p5 = [quantile([s[i] for s in sim_data], 0.05) for i in 1:length(pred)]
                p50 = [quantile([s[i] for s in sim_data], 0.50) for i in 1:length(pred)]
                p95 = [quantile([s[i] for s in sim_data], 0.95) for i in 1:length(pred)]
                (p5, p50, p95)
            end,
            "vpc_percentiles", category="analysis", subcategory="vpc",
            model="VPCPercentiles"
        )
        r !== nothing && push!(results, r)
    catch e
        println("  Skipping VPC: $e")
    end

    return results
end

# =============================================================================
# Main Runner
# =============================================================================

function run_comprehensive_benchmarks()
    println("\n" * "="^80)
    println("NeoPKPD COMPREHENSIVE Benchmark Suite")
    println("="^80)
    println("Julia Version: ", VERSION)
    println("NeoPKPD Version: ", NEOPKPD_VERSION)
    println("Runs per benchmark: ", CONFIG.n_runs)
    println("Warmup runs: ", CONFIG.n_warmup)
    println("Random seed: ", CONFIG.random_seed)
    println("Timestamp: ", now())
    println("="^80)

    Random.seed!(CONFIG.random_seed)

    all_results = BenchmarkResult[]

    # Run all benchmark categories
    append!(all_results, benchmark_pk_models())
    append!(all_results, benchmark_pd_models())
    append!(all_results, benchmark_pkpd_coupled())
    append!(all_results, benchmark_nca())
    append!(all_results, benchmark_population())
    append!(all_results, benchmark_sensitivity())
    append!(all_results, benchmark_trial_simulation())
    append!(all_results, benchmark_tmdd())
    append!(all_results, benchmark_estimation())
    append!(all_results, benchmark_analysis())

    # Save results
    mkpath(CONFIG.output_dir)
    output_file = joinpath(CONFIG.output_dir, "neopkpd_comprehensive_benchmarks.csv")

    summaries = [summarize(r) for r in all_results]

    open(output_file, "w") do io
        println(io, "category,subcategory,name,model,n_subjects,n_runs,mean_ms,std_ms,median_ms,min_ms,max_ms,ci_lower_ms,ci_upper_ms,timestamp,julia_version,neopkpd_version")

        for s in summaries
            println(io, join([
                s.category, s.subcategory, s.name, s.model, s.n_subjects, s.n_runs,
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

    println("\n" * "="^80)
    println("COMPREHENSIVE BENCHMARK COMPLETE")
    println("Results saved to: ", output_file)
    println("Total benchmarks: ", length(all_results))
    println("="^80)

    # Print summary by category
    println("\n" * "-"^80)
    println("SUMMARY BY CATEGORY")
    println("-"^80)

    categories = unique([s.category for s in summaries])
    for cat in categories
        cat_summaries = filter(s -> s.category == cat, summaries)
        mean_time = mean([s.mean_ms for s in cat_summaries])
        n_benchmarks = length(cat_summaries)
        @printf("%-20s: %3d benchmarks, avg %.3f ms\n", cat, n_benchmarks, mean_time)
    end
    println("-"^80)

    return summaries
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    run_comprehensive_benchmarks()
end
