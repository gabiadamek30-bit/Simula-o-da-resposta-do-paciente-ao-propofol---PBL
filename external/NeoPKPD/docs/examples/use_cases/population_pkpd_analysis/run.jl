# Population PKPD Analysis - Complete Workflow
# Run: julia --project=packages/core docs/examples/use_cases/population_pkpd_analysis/run.jl

using Pkg
Pkg.activate("packages/core")

using NeoPKPD
using JSON
using Statistics
using Random

const BASE_DIR = "docs/examples/use_cases/population_pkpd_analysis"
const N_SUBJECTS = 50
const DOSE = 50.0  # mg

# ============================================================================
# Step 1: Define Model Parameters
# ============================================================================

# PK Parameters (Two-compartment oral)
const PK_PARAMS = Dict(
    :Ka => 1.0,   # /h
    :CL => 5.0,   # L/h
    :V1 => 20.0,  # L
    :Q => 2.0,    # L/h
    :V2 => 40.0   # L
)

# PD Parameters (Indirect response - inhibition of production)
const PD_PARAMS = Dict(
    :Kin => 1.0,   # /h
    :Kout => 0.1,  # /h
    :Imax => 0.9,  # Maximal inhibition
    :IC50 => 2.0   # mg/L
)

# IIV (CV%)
const OMEGA = Dict(
    :Ka => 0.30,
    :CL => 0.25,
    :V1 => 0.20,
    :Kin => 0.20,
    :Kout => 0.20,
    :Imax => 0.15,
    :IC50 => 0.30
)

# ============================================================================
# Step 2: Generate Population
# ============================================================================

function generate_population(seed::UInt64)
    Random.seed!(seed)

    subjects = Vector{Dict{String, Any}}()

    for i in 1:N_SUBJECTS
        # Weight covariate (50-100 kg)
        wt = 50.0 + 50.0 * rand()

        # Individual random effects
        etas = Dict{Symbol, Float64}()
        for (param, omega) in OMEGA
            etas[param] = randn() * omega
        end

        # Individual PK parameters with weight scaling
        cl_ind = PK_PARAMS[:CL] * (wt/70)^0.75 * exp(etas[:CL])
        v1_ind = PK_PARAMS[:V1] * (wt/70)^1.0 * exp(etas[:V1])
        ka_ind = PK_PARAMS[:Ka] * exp(etas[:Ka])

        # Individual PD parameters
        kin_ind = PD_PARAMS[:Kin] * exp(etas[:Kin])
        kout_ind = PD_PARAMS[:Kout] * exp(etas[:Kout])
        imax_ind = PD_PARAMS[:Imax] * exp(etas[:Imax])
        ic50_ind = PD_PARAMS[:IC50] * exp(etas[:IC50])

        push!(subjects, Dict(
            "id" => i,
            "wt" => wt,
            "pk" => Dict(
                "Ka" => ka_ind,
                "CL" => cl_ind,
                "V1" => v1_ind,
                "Q" => PK_PARAMS[:Q],
                "V2" => PK_PARAMS[:V2]
            ),
            "pd" => Dict(
                "Kin" => kin_ind,
                "Kout" => kout_ind,
                "Imax" => imax_ind,
                "IC50" => ic50_ind,
                "E0" => kin_ind / kout_ind  # Baseline
            ),
            "etas" => etas
        ))
    end

    return subjects
end

# ============================================================================
# Step 3: Simulate PK
# ============================================================================

function simulate_pk(subjects::Vector)
    println("  Simulating PK profiles...")

    pk_grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pk_results = Vector{Dict{String, Any}}()

    for subj in subjects
        p = subj["pk"]

        model = ModelSpec(
            TwoCompIVBolus(),  # Using IV for simplicity; oral would need absorption
            "pkpd_pk_$(subj["id"])",
            TwoCompIVBolusParams(p["CL"], p["V1"], p["Q"], p["V2"]),
            [DoseEvent(0.0, DOSE)]
        )

        res = simulate(model, pk_grid, solver)

        push!(pk_results, Dict(
            "id" => subj["id"],
            "times" => collect(pk_grid.saveat),
            "conc" => res.observations[:conc],
            "params" => p
        ))
    end

    return pk_results
end

# ============================================================================
# Step 4: Simulate PD (Indirect Response)
# ============================================================================

function simulate_pd_indirect(subjects::Vector, pk_results::Vector)
    println("  Simulating PD (indirect response)...")

    # PD needs longer simulation (72h) to see full response
    pd_times = collect(0.0:1.0:72.0)

    pd_results = Vector{Dict{String, Any}}()

    for (subj, pk_res) in zip(subjects, pk_results)
        pd_params = subj["pd"]

        # Interpolate PK concentrations for PD simulation
        pk_times = pk_res["times"]
        pk_conc = pk_res["conc"]

        # Extend PK with elimination phase
        function get_conc(t)
            if t <= pk_times[end]
                # Linear interpolation
                idx = searchsortedlast(pk_times, t)
                if idx == 0
                    return pk_conc[1]
                elseif idx >= length(pk_times)
                    return pk_conc[end]
                else
                    frac = (t - pk_times[idx]) / (pk_times[idx+1] - pk_times[idx])
                    return pk_conc[idx] + frac * (pk_conc[idx+1] - pk_conc[idx])
                end
            else
                # Extrapolate with terminal elimination
                return pk_conc[end] * exp(-0.1 * (t - pk_times[end]))
            end
        end

        # Simulate indirect response
        # dE/dt = Kin * (1 - Imax * C / (IC50 + C)) - Kout * E
        E = zeros(length(pd_times))
        E[1] = pd_params["E0"]  # Baseline

        for i in 2:length(pd_times)
            dt = pd_times[i] - pd_times[i-1]
            t = pd_times[i-1]
            C = max(get_conc(t), 0.0)

            # Inhibition term
            inhibition = pd_params["Imax"] * C / (pd_params["IC50"] + C)

            # Euler integration
            dE = pd_params["Kin"] * (1 - inhibition) - pd_params["Kout"] * E[i-1]
            E[i] = max(E[i-1] + dE * dt, 0.0)
        end

        push!(pd_results, Dict(
            "id" => subj["id"],
            "times" => pd_times,
            "effect" => E,
            "baseline" => pd_params["E0"],
            "params" => pd_params
        ))
    end

    return pd_results
end

# ============================================================================
# Step 5: Compute Summary Statistics
# ============================================================================

function compute_summary(pk_results::Vector, pd_results::Vector)
    println("  Computing summary statistics...")

    # PK metrics
    pk_cmax = [maximum(r["conc"]) for r in pk_results]
    pk_tmax = [r["times"][argmax(r["conc"])] for r in pk_results]

    # Compute AUC for each subject
    pk_auc = Float64[]
    for r in pk_results
        auc = 0.0
        t = r["times"]
        c = r["conc"]
        for i in 2:length(t)
            auc += 0.5 * (c[i] + c[i-1]) * (t[i] - t[i-1])
        end
        push!(pk_auc, auc)
    end

    # PD metrics
    pd_min = [minimum(r["effect"]) for r in pd_results]
    pd_baseline = [r["baseline"] for r in pd_results]
    pd_max_inhibition = [(1 - minimum(r["effect"]) / r["baseline"]) * 100 for r in pd_results]
    pd_tmin = [r["times"][argmin(r["effect"])] for r in pd_results]

    summary = Dict(
        "n" => length(pk_results),
        "pk" => Dict(
            "Cmax_mean" => mean(pk_cmax),
            "Cmax_sd" => std(pk_cmax),
            "Tmax_median" => median(pk_tmax),
            "AUC_mean" => mean(pk_auc),
            "AUC_sd" => std(pk_auc)
        ),
        "pd" => Dict(
            "baseline_mean" => mean(pd_baseline),
            "baseline_sd" => std(pd_baseline),
            "max_inhibition_mean" => mean(pd_max_inhibition),
            "max_inhibition_sd" => std(pd_max_inhibition),
            "time_to_nadir_mean" => mean(pd_tmin),
            "time_to_nadir_sd" => std(pd_tmin)
        )
    )

    return summary
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("=" ^ 60)
    println("Population PKPD Analysis")
    println("=" ^ 60)

    # Step 1: Generate population
    println("\nStep 1: Generating $(N_SUBJECTS) subjects...")
    subjects = generate_population(UInt64(20240201))

    weights = [s["wt"] for s in subjects]
    println("  Weight range: $(round(minimum(weights), digits=1)) - $(round(maximum(weights), digits=1)) kg")

    # Step 2: Simulate PK
    println("\nStep 2: Running simulations...")
    pk_results = simulate_pk(subjects)

    # Step 3: Simulate PD
    pd_results = simulate_pd_indirect(subjects, pk_results)

    # Step 4: Compute summary
    println("\nStep 3: Computing summary...")
    summary = compute_summary(pk_results, pd_results)

    # Print summary
    println("\n" * "-" ^ 40)
    println("POPULATION PKPD SUMMARY")
    println("-" ^ 40)
    println("N = $(summary["n"])")
    println()
    println("PK Metrics:")
    println("  Cmax: $(round(summary["pk"]["Cmax_mean"], digits=2)) +/- $(round(summary["pk"]["Cmax_sd"], digits=2)) mg/L")
    println("  Tmax (median): $(round(summary["pk"]["Tmax_median"], digits=1)) h")
    println("  AUC_0_24: $(round(summary["pk"]["AUC_mean"], digits=1)) +/- $(round(summary["pk"]["AUC_sd"], digits=1)) mg*h/L")
    println()
    println("PD Metrics:")
    println("  Baseline E0: $(round(summary["pd"]["baseline_mean"], digits=2)) +/- $(round(summary["pd"]["baseline_sd"], digits=2))")
    println("  Max inhibition: $(round(summary["pd"]["max_inhibition_mean"], digits=1)) +/- $(round(summary["pd"]["max_inhibition_sd"], digits=1)) %")
    println("  Time to nadir: $(round(summary["pd"]["time_to_nadir_mean"], digits=1)) +/- $(round(summary["pd"]["time_to_nadir_sd"], digits=1)) h")
    println("-" ^ 40)

    # Step 5: Save outputs
    println("\nStep 4: Saving outputs...")
    mkpath(joinpath(BASE_DIR, "output"))

    # PK results
    pk_path = joinpath(BASE_DIR, "output", "pk_results.json")
    open(pk_path, "w") do io
        JSON.print(io, pk_results, 2)
    end
    println("  Wrote: $pk_path")

    # PD results
    pd_path = joinpath(BASE_DIR, "output", "pd_results.json")
    open(pd_path, "w") do io
        JSON.print(io, pd_results, 2)
    end
    println("  Wrote: $pd_path")

    # Summary
    summary_path = joinpath(BASE_DIR, "output", "summary.json")
    open(summary_path, "w") do io
        JSON.print(io, summary, 2)
    end
    println("  Wrote: $summary_path")

    println("\n" * "=" ^ 60)
    println("Population PKPD analysis complete!")
    println("=" ^ 60)
end

main()
