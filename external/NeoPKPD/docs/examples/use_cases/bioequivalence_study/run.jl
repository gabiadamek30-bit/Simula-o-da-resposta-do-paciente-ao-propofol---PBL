# Bioequivalence Study Simulation and Analysis
# Run: julia --project=packages/core docs/examples/use_cases/bioequivalence_study/run.jl

using Pkg
Pkg.activate("packages/core")

using NeoPKPD
using JSON
using Statistics
using Random

const BASE_DIR = "docs/examples/use_cases/bioequivalence_study"

# Study design parameters
const N_SUBJECTS_PER_SEQUENCE = 12
const N_TOTAL = 2 * N_SUBJECTS_PER_SEQUENCE

# ============================================================================
# Step 1: Define Model Parameters
# ============================================================================

struct FormulationParams
    name::String
    ka::Float64
    cl::Float64
    v::Float64
end

const REFERENCE = FormulationParams("Reference", 1.5, 10.0, 100.0)
const TEST = FormulationParams("Test", 1.4, 10.0, 100.0)  # Slightly slower absorption

# IIV parameters (same for both formulations - within-subject design)
const OMEGA_KA = 0.30  # 30% CV
const OMEGA_CL = 0.25  # 25% CV
const OMEGA_V = 0.20   # 20% CV

# Dose
const DOSE = 100.0  # mg

# ============================================================================
# Step 2: Generate Population with Consistent Random Effects
# ============================================================================

function generate_population(seed::UInt64)
    Random.seed!(seed)

    subjects = Vector{Dict{String, Any}}()

    for i in 1:N_TOTAL
        # Generate individual random effects (eta)
        eta_ka = randn() * OMEGA_KA
        eta_cl = randn() * OMEGA_CL
        eta_v = randn() * OMEGA_V

        # Sequence assignment (1 = TR, 2 = RT)
        sequence = i <= N_SUBJECTS_PER_SEQUENCE ? 1 : 2

        push!(subjects, Dict(
            "id" => i,
            "sequence" => sequence,
            "eta_ka" => eta_ka,
            "eta_cl" => eta_cl,
            "eta_v" => eta_v
        ))
    end

    return subjects
end

# ============================================================================
# Step 3: Simulate One Period
# ============================================================================

function simulate_period(subjects::Vector, formulation::FormulationParams, period::Int)
    println("  Simulating Period $period: $(formulation.name)")

    results = Vector{Dict{String, Any}}()

    grid = SimGrid(0.0, 24.0, [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0])
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    for subj in subjects
        # Apply individual random effects
        ka_ind = formulation.ka * exp(subj["eta_ka"])
        cl_ind = formulation.cl * exp(subj["eta_cl"])
        v_ind = formulation.v * exp(subj["eta_v"])

        # Create model spec
        model = ModelSpec(
            OneCompOralFirstOrder(),
            "be_$(formulation.name)_$(subj["id"])",
            OneCompOralFirstOrderParams(ka_ind, cl_ind, v_ind),
            [DoseEvent(0.0, DOSE)]
        )

        # Simulate
        res = simulate(model, grid, solver)

        # Extract concentrations
        conc = res.observations[:conc]

        push!(results, Dict(
            "id" => subj["id"],
            "sequence" => subj["sequence"],
            "period" => period,
            "formulation" => formulation.name,
            "times" => collect(grid.saveat),
            "conc" => conc,
            "params" => Dict("Ka" => ka_ind, "CL" => cl_ind, "V" => v_ind)
        ))
    end

    return results
end

# ============================================================================
# Step 4: NCA Analysis
# ============================================================================

function compute_nca(result::Dict)
    times = result["times"]
    conc = result["conc"]

    # Cmax and Tmax
    cmax_idx = argmax(conc)
    cmax = conc[cmax_idx]
    tmax = times[cmax_idx]

    # AUC by trapezoidal rule
    auc = 0.0
    for i in 2:length(times)
        dt = times[i] - times[i-1]
        auc += 0.5 * (conc[i] + conc[i-1]) * dt
    end

    return Dict(
        "Cmax" => cmax,
        "Tmax" => tmax,
        "AUC_0_24" => auc
    )
end

# ============================================================================
# Step 5: Statistical Analysis
# ============================================================================

function compute_be_statistics(all_results::Vector)
    # Organize by subject
    subject_data = Dict{Int, Dict{String, Any}}()

    for res in all_results
        id = res["id"]
        if !haskey(subject_data, id)
            subject_data[id] = Dict{String, Any}("sequence" => res["sequence"])
        end

        nca = compute_nca(res)
        form = res["formulation"]

        subject_data[id]["$(form)_Cmax"] = nca["Cmax"]
        subject_data[id]["$(form)_AUC"] = nca["AUC_0_24"]
    end

    # Compute log ratios (T/R) for each subject
    log_cmax_ratios = Float64[]
    log_auc_ratios = Float64[]

    for (id, data) in subject_data
        if haskey(data, "Test_Cmax") && haskey(data, "Reference_Cmax")
            push!(log_cmax_ratios, log(data["Test_Cmax"]) - log(data["Reference_Cmax"]))
            push!(log_auc_ratios, log(data["Test_AUC"]) - log(data["Reference_AUC"]))
        end
    end

    n = length(log_cmax_ratios)

    # Geometric mean ratios
    gmr_cmax = exp(mean(log_cmax_ratios))
    gmr_auc = exp(mean(log_auc_ratios))

    # 90% CI (using t-distribution approximation)
    t_crit = 1.714  # t(0.95, 22) for n=24 in 2x2 crossover

    se_cmax = std(log_cmax_ratios) / sqrt(n)
    ci_cmax_lower = exp(mean(log_cmax_ratios) - t_crit * se_cmax)
    ci_cmax_upper = exp(mean(log_cmax_ratios) + t_crit * se_cmax)

    se_auc = std(log_auc_ratios) / sqrt(n)
    ci_auc_lower = exp(mean(log_auc_ratios) - t_crit * se_auc)
    ci_auc_upper = exp(mean(log_auc_ratios) + t_crit * se_auc)

    # BE assessment
    cmax_be = ci_cmax_lower >= 0.80 && ci_cmax_upper <= 1.25
    auc_be = ci_auc_lower >= 0.80 && ci_auc_upper <= 1.25

    return Dict(
        "n_subjects" => n,
        "Cmax" => Dict(
            "GMR" => gmr_cmax,
            "CI_lower" => ci_cmax_lower,
            "CI_upper" => ci_cmax_upper,
            "BE" => cmax_be
        ),
        "AUC" => Dict(
            "GMR" => gmr_auc,
            "CI_lower" => ci_auc_lower,
            "CI_upper" => ci_auc_upper,
            "BE" => auc_be
        ),
        "overall_BE" => cmax_be && auc_be
    )
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("=" ^ 60)
    println("Bioequivalence Study Simulation")
    println("=" ^ 60)

    # Step 1: Generate population
    println("\nStep 1: Generating $(N_TOTAL) subjects...")
    subjects = generate_population(UInt64(20240101))
    println("  Sequence 1 (T-R): $(N_SUBJECTS_PER_SEQUENCE) subjects")
    println("  Sequence 2 (R-T): $(N_SUBJECTS_PER_SEQUENCE) subjects")

    # Step 2: Simulate crossover study
    println("\nStep 2: Simulating crossover study...")
    all_results = Vector{Dict{String, Any}}()

    # Period 1
    # Sequence 1 gets Test, Sequence 2 gets Reference
    seq1_period1 = [s for s in subjects if s["sequence"] == 1]
    seq2_period1 = [s for s in subjects if s["sequence"] == 2]

    append!(all_results, simulate_period(seq1_period1, TEST, 1))
    append!(all_results, simulate_period(seq2_period1, REFERENCE, 1))

    # Period 2 (crossover)
    # Sequence 1 gets Reference, Sequence 2 gets Test
    append!(all_results, simulate_period(seq1_period1, REFERENCE, 2))
    append!(all_results, simulate_period(seq2_period1, TEST, 2))

    # Step 3: NCA and statistical analysis
    println("\nStep 3: Computing NCA metrics and BE statistics...")
    be_stats = compute_be_statistics(all_results)

    # Print results
    println("\n" * "-" ^ 40)
    println("BIOEQUIVALENCE RESULTS")
    println("-" ^ 40)
    println("Subjects evaluated: $(be_stats["n_subjects"])")
    println()
    println("Cmax:")
    println("  GMR (T/R): $(round(be_stats["Cmax"]["GMR"], digits=3))")
    println("  90% CI: $(round(be_stats["Cmax"]["CI_lower"], digits=3)) - $(round(be_stats["Cmax"]["CI_upper"], digits=3))")
    println("  BE: $(be_stats["Cmax"]["BE"] ? "PASS" : "FAIL")")
    println()
    println("AUC_0_24:")
    println("  GMR (T/R): $(round(be_stats["AUC"]["GMR"], digits=3))")
    println("  90% CI: $(round(be_stats["AUC"]["CI_lower"], digits=3)) - $(round(be_stats["AUC"]["CI_upper"], digits=3))")
    println("  BE: $(be_stats["AUC"]["BE"] ? "PASS" : "FAIL")")
    println()
    println("Overall BE Conclusion: $(be_stats["overall_BE"] ? "BIOEQUIVALENT" : "NOT BIOEQUIVALENT")")
    println("-" ^ 40)

    # Step 4: Save outputs
    println("\nStep 4: Saving outputs...")
    mkpath(joinpath(BASE_DIR, "output"))

    # Individual results
    ind_path = joinpath(BASE_DIR, "output", "individual_results.json")
    open(ind_path, "w") do io
        JSON.print(io, all_results, 2)
    end
    println("  Wrote: $ind_path")

    # BE statistics
    stats_path = joinpath(BASE_DIR, "output", "be_statistics.json")
    open(stats_path, "w") do io
        JSON.print(io, be_stats, 2)
    end
    println("  Wrote: $stats_path")

    println("\n" * "=" ^ 60)
    println("Bioequivalence study complete!")
    println("=" ^ 60)
end

main()
