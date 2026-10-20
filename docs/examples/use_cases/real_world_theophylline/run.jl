# Real-World Theophylline Analysis - Complete Workflow
# Run: julia --project=packages/core docs/examples/use_cases/real_world_theophylline/run.jl

using Pkg
Pkg.activate("packages/core")

using NeoPKPD
using JSON
using Statistics
using DelimitedFiles

const BASE_DIR = "docs/examples/use_cases/real_world_theophylline"
const DATA_FILE = "docs/examples/real_world_validation/datasets/theophylline_theo_sd/theo_sd.csv"

# ============================================================================
# Step 1: Load and Parse Data
# ============================================================================

function load_theophylline_data()
    println("Step 1: Loading theophylline data...")

    data, header = readdlm(DATA_FILE, ',', Any, header=true)
    header = vec(header)

    col_idx = Dict(String(h) => i for (i, h) in enumerate(header))

    subjects = Dict{Int, Dict{String, Any}}()

    for row in eachrow(data)
        id = Int(row[col_idx["ID"]])
        time = Float64(row[col_idx["TIME"]])
        dv = Float64(row[col_idx["DV"]])
        amt = Float64(row[col_idx["AMT"]])
        evid = Int(row[col_idx["EVID"]])
        wt = Float64(row[col_idx["WT"]])

        if !haskey(subjects, id)
            subjects[id] = Dict(
                "id" => id,
                "wt" => wt,
                "dose" => 0.0,
                "times" => Float64[],
                "conc" => Float64[]
            )
        end

        if evid == 101  # Dose event
            subjects[id]["dose"] = amt
        elseif evid == 0 && time > 0  # Observation (exclude baseline)
            push!(subjects[id]["times"], time)
            push!(subjects[id]["conc"], dv)
        end
    end

    println("  Loaded $(length(subjects)) subjects")
    return subjects
end

# ============================================================================
# Step 2: NCA Analysis
# ============================================================================

function compute_nca_metrics(times::Vector{Float64}, conc::Vector{Float64}, dose::Float64)
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

    # Terminal half-life (using last 3 points)
    if length(times) >= 3
        n_term = min(4, length(times) - cmax_idx)
        if n_term >= 2
            term_idx = (length(times) - n_term + 1):length(times)
            t_term = times[term_idx]
            c_term = log.(max.(conc[term_idx], 1e-10))

            # Linear regression for log(C) vs t
            n = length(t_term)
            sum_t = sum(t_term)
            sum_c = sum(c_term)
            sum_tc = sum(t_term .* c_term)
            sum_t2 = sum(t_term .^ 2)

            lambda_z = -(n * sum_tc - sum_t * sum_c) / (n * sum_t2 - sum_t^2)
            t_half = log(2) / max(lambda_z, 1e-10)
        else
            t_half = NaN
        end
    else
        t_half = NaN
    end

    return Dict(
        "Cmax" => cmax,
        "Tmax" => tmax,
        "AUC_0_t" => auc,
        "t_half" => t_half,
        "dose" => dose
    )
end

function run_nca_analysis(subjects::Dict)
    println("\nStep 2: Running NCA analysis...")

    nca_results = Dict{Int, Dict}()

    for (id, subj) in subjects
        nca = compute_nca_metrics(subj["times"], subj["conc"], subj["dose"])
        nca["id"] = id
        nca["wt"] = subj["wt"]
        nca_results[id] = nca

        println("  Subject $id: Cmax=$(round(nca["Cmax"], digits=2)), Tmax=$(round(nca["Tmax"], digits=2)), AUC=$(round(nca["AUC_0_t"], digits=1))")
    end

    # Population summary
    cmax_all = [nca["Cmax"] for nca in values(nca_results)]
    tmax_all = [nca["Tmax"] for nca in values(nca_results)]
    auc_all = [nca["AUC_0_t"] for nca in values(nca_results)]

    summary = Dict(
        "n" => length(nca_results),
        "Cmax_mean" => mean(cmax_all),
        "Cmax_sd" => std(cmax_all),
        "Tmax_median" => median(tmax_all),
        "AUC_mean" => mean(auc_all),
        "AUC_sd" => std(auc_all)
    )

    println("\n  Population Summary:")
    println("    Cmax: $(round(summary["Cmax_mean"], digits=2)) +/- $(round(summary["Cmax_sd"], digits=2)) mg/L")
    println("    Tmax (median): $(round(summary["Tmax_median"], digits=2)) h")
    println("    AUC: $(round(summary["AUC_mean"], digits=1)) +/- $(round(summary["AUC_sd"], digits=1)) mg*h/L")

    return nca_results, summary
end

# ============================================================================
# Step 3: Population Simulation with Fitted Parameters
# ============================================================================

function run_population_simulation(subjects::Dict)
    println("\nStep 3: Running population simulation...")

    # Typical parameters from literature for theophylline oral
    # Ka ~1.5 /h, CL ~2.8 L/h (0.04 L/h/kg * 70kg), V ~35 L (0.5 L/kg * 70kg)
    typical_ka = 1.5
    typical_cl = 2.8  # L/h
    typical_v = 35.0  # L

    # Mean dose and weight
    doses = [subj["dose"] for subj in values(subjects)]
    weights = [subj["wt"] for subj in values(subjects)]
    mean_dose = mean(doses)

    # Create model spec
    model = ModelSpec(
        OneCompOralFirstOrder(),
        "theophylline_oral",
        OneCompOralFirstOrderParams(typical_ka, typical_cl, typical_v),
        [DoseEvent(0.0, mean_dose)]
    )

    # Setup covariate model (allometric scaling)
    cov_model = CovariateModel(
        "wt_scaling",
        [
            CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),
            CovariateEffect(PowerCovariate(), :V, :WT, 1.0, 70.0),
        ]
    )

    # Individual covariates from actual subjects
    ind_covs = IndividualCovariates[]
    for id in sort(collect(keys(subjects)))
        subj = subjects[id]
        push!(ind_covs, IndividualCovariates(Dict(:WT => subj["wt"]), nothing))
    end

    # IIV specification
    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.20, :V => 0.15, :Ka => 0.25), UInt64(12345), length(ind_covs))

    # Population spec
    pop = PopulationSpec(model, iiv, nothing, cov_model, ind_covs)

    # Simulation grid
    grid = SimGrid(0.0, 25.0, collect(0.0:0.25:25.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    # Run simulation
    result = simulate_population(pop, grid, solver)

    println("  Simulated $(length(result.individuals)) subjects")
    println("  Time points: $(length(grid.saveat))")

    return pop, grid, solver, result
end

# ============================================================================
# Step 4: VPC Calculation
# ============================================================================

function compute_vpc_statistics(subjects::Dict, pop_result)
    println("\nStep 4: Computing VPC statistics...")

    # Bin observed data
    time_bins = [0.0, 1.0, 2.0, 4.0, 6.0, 9.0, 12.0, 24.0]

    # Collect observed concentrations by bin
    obs_by_bin = Dict{Int, Vector{Float64}}()
    for i in 1:(length(time_bins)-1)
        obs_by_bin[i] = Float64[]
    end

    for subj in values(subjects)
        for (t, c) in zip(subj["times"], subj["conc"])
            for i in 1:(length(time_bins)-1)
                if time_bins[i] <= t < time_bins[i+1]
                    push!(obs_by_bin[i], c)
                    break
                end
            end
        end
    end

    # Compute observed percentiles
    obs_percentiles = Dict{Int, Dict{String, Float64}}()
    for (bin, concs) in obs_by_bin
        if length(concs) >= 2
            obs_percentiles[bin] = Dict(
                "p5" => quantile(concs, 0.05),
                "p50" => quantile(concs, 0.50),
                "p95" => quantile(concs, 0.95),
                "n" => Float64(length(concs))
            )
        end
    end

    # Compute simulated percentiles from population result
    sim_conc = pop_result.summaries[:conc]

    # Extract quantiles from the Dict
    sim_p5 = get(sim_conc.quantiles, 0.05, sim_conc.mean)
    sim_p95 = get(sim_conc.quantiles, 0.95, sim_conc.mean)

    vpc_summary = Dict(
        "time_bins" => time_bins,
        "observed" => obs_percentiles,
        "simulated_mean" => sim_conc.mean,
        "simulated_p5" => sim_p5,
        "simulated_p95" => sim_p95
    )

    println("  Computed VPC statistics for $(length(obs_percentiles)) time bins")

    return vpc_summary
end

# ============================================================================
# Step 5: Generate Outputs
# ============================================================================

function save_outputs(nca_results, nca_summary, pop, grid, solver, pop_result, vpc_summary)
    println("\nStep 5: Saving outputs...")

    mkpath(joinpath(BASE_DIR, "output"))

    # NCA results
    nca_output = Dict(
        "individual" => [nca_results[id] for id in sort(collect(keys(nca_results)))],
        "summary" => nca_summary
    )
    nca_path = joinpath(BASE_DIR, "output", "nca_results.json")
    open(nca_path, "w") do io
        JSON.print(io, nca_output, 2)
    end
    println("  Wrote: $nca_path")

    # Population artifact
    pop_path = joinpath(BASE_DIR, "output", "population_simulation.json")
    write_population_json(pop_path; population_spec=pop, grid=grid, solver=solver, result=pop_result)
    println("  Wrote: $pop_path")

    # VPC summary
    vpc_path = joinpath(BASE_DIR, "output", "vpc_summary.json")
    open(vpc_path, "w") do io
        JSON.print(io, vpc_summary, 2)
    end
    println("  Wrote: $vpc_path")

    # Summary report
    report = Dict(
        "analysis" => "Theophylline Single-Dose PK Analysis",
        "dataset" => "theo_sd (nlmixr2)",
        "n_subjects" => nca_summary["n"],
        "nca_summary" => Dict(
            "Cmax_mean_mg_L" => round(nca_summary["Cmax_mean"], digits=2),
            "Cmax_sd" => round(nca_summary["Cmax_sd"], digits=2),
            "Tmax_median_h" => round(nca_summary["Tmax_median"], digits=2),
            "AUC_mean_mg_h_L" => round(nca_summary["AUC_mean"], digits=1),
            "AUC_sd" => round(nca_summary["AUC_sd"], digits=1)
        ),
        "model" => "OneCompOralFirstOrder with weight-based allometric scaling",
        "parameters" => Dict(
            "Ka" => "1.5 /h",
            "CL" => "2.8 L/h (0.04 L/h/kg)",
            "V" => "35 L (0.5 L/kg)"
        )
    )
    report_path = joinpath(BASE_DIR, "output", "report.json")
    open(report_path, "w") do io
        JSON.print(io, report, 2)
    end
    println("  Wrote: $report_path")
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("=" ^ 60)
    println("Real-World Theophylline PK Analysis")
    println("=" ^ 60)

    # Step 1: Load data
    subjects = load_theophylline_data()

    # Step 2: NCA analysis
    nca_results, nca_summary = run_nca_analysis(subjects)

    # Step 3: Population simulation
    pop, grid, solver, pop_result = run_population_simulation(subjects)

    # Step 4: VPC statistics
    vpc_summary = compute_vpc_statistics(subjects, pop_result)

    # Step 5: Save outputs
    save_outputs(nca_results, nca_summary, pop, grid, solver, pop_result, vpc_summary)

    println("\n" * "=" ^ 60)
    println("Analysis complete!")
    println("=" ^ 60)
end

main()
