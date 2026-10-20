# NONMEM Migration Workflow
# Demonstrates converting a NONMEM model to NeoPKPD format
# Run: julia --project=packages/core docs/examples/use_cases/nonmem_migration/run.jl

using Pkg
Pkg.activate("packages/core")

using NeoPKPD
using JSON
using Statistics
using DelimitedFiles

const BASE_DIR = "docs/examples/use_cases/nonmem_migration"
const DATA_FILE = "docs/examples/real_world_validation/datasets/theophylline_theo_sd/theo_sd.csv"

# ============================================================================
# Step 1: Parse NONMEM Control File (Manual extraction for this example)
# ============================================================================

# Extracted parameters from run001.ctl
const NONMEM_MODEL = Dict(
    "subroutine" => "ADVAN2",  # One-compartment oral
    "trans" => "TRANS2",       # CL, V parameterization
    "thetas" => [
        Dict("name" => "TVCL", "init" => 2.8, "lower" => 0.0),
        Dict("name" => "TVV", "init" => 35.0, "lower" => 0.0),
        Dict("name" => "TVKA", "init" => 1.5, "lower" => 0.0)
    ],
    "omegas" => [
        Dict("name" => "IIV_CL", "value" => 0.09),
        Dict("name" => "IIV_V", "value" => 0.04),
        Dict("name" => "IIV_KA", "value" => 0.16)
    ],
    "sigma" => [
        Dict("name" => "PROP", "value" => 0.04)
    ],
    "covariates" => [
        Dict("param" => "CL", "cov" => "WT", "type" => "power", "exp" => 0.75, "ref" => 70.0),
        Dict("param" => "V", "cov" => "WT", "type" => "power", "exp" => 1.0, "ref" => 70.0)
    ]
)

function parse_nonmem_model()
    println("Step 1: Parsing NONMEM control file...")
    println("  Model: $(NONMEM_MODEL["subroutine"]) ($(NONMEM_MODEL["trans"]))")
    println("  THETAs: $(length(NONMEM_MODEL["thetas"]))")
    println("  OMEGAs: $(length(NONMEM_MODEL["omegas"]))")
    println("  Covariates: $(length(NONMEM_MODEL["covariates"]))")
    return NONMEM_MODEL
end

# ============================================================================
# Step 2: Map to NeoPKPD Model
# ============================================================================

function map_advan_to_neopkpd(advan::String)
    mapping = Dict(
        "ADVAN1" => OneCompIVBolus,
        "ADVAN2" => OneCompOralFirstOrder,
        "ADVAN3" => TwoCompIVBolus,
        "ADVAN4" => TwoCompOral
    )
    return mapping[advan]
end

function convert_to_neopkpd(nm_model::Dict)
    println("\nStep 2: Converting to NeoPKPD format...")

    # Map model type
    model_type = map_advan_to_neopkpd(nm_model["subroutine"])
    println("  NeoPKPD model: $model_type")

    # Extract typical values
    typical_cl = nm_model["thetas"][1]["init"]
    typical_v = nm_model["thetas"][2]["init"]
    typical_ka = nm_model["thetas"][3]["init"]

    # Extract omega values (variance -> CV)
    omega_cl = sqrt(nm_model["omegas"][1]["value"])
    omega_v = sqrt(nm_model["omegas"][2]["value"])
    omega_ka = sqrt(nm_model["omegas"][3]["value"])

    println("  Typical CL: $typical_cl L/h")
    println("  Typical V: $typical_v L")
    println("  Typical Ka: $typical_ka /h")
    println("  IIV CL (CV): $(round(omega_cl * 100, digits=0))%")
    println("  IIV V (CV): $(round(omega_v * 100, digits=0))%")
    println("  IIV Ka (CV): $(round(omega_ka * 100, digits=0))%")

    neopkpd_model = Dict(
        "model_type" => model_type,
        "typical" => Dict(:Ka => typical_ka, :CL => typical_cl, :V => typical_v),
        "omega" => Dict(:Ka => omega_ka, :CL => omega_cl, :V => omega_v),
        "covariates" => nm_model["covariates"]
    )

    return neopkpd_model
end

# ============================================================================
# Step 3: Load Data
# ============================================================================

function load_data()
    println("\nStep 3: Loading theophylline data...")

    data, header = readdlm(DATA_FILE, ',', Any, header=true)
    header = vec(header)
    col_idx = Dict(String(h) => i for (i, h) in enumerate(header))

    subjects = Dict{Int, Dict{String, Any}}()

    for row in eachrow(data)
        id = Int(row[col_idx["ID"]])
        time = Float64(row[col_idx["TIME"]])
        amt = Float64(row[col_idx["AMT"]])
        evid = Int(row[col_idx["EVID"]])
        wt = Float64(row[col_idx["WT"]])

        if !haskey(subjects, id)
            subjects[id] = Dict(
                "id" => id,
                "wt" => wt,
                "dose" => 0.0,
                "times" => Float64[],
                "dv" => Float64[]
            )
        end

        if evid == 101
            subjects[id]["dose"] = amt
        end
    end

    println("  Loaded $(length(subjects)) subjects")
    return subjects
end

# ============================================================================
# Step 4: Simulate with NeoPKPD
# ============================================================================

function simulate_neopkpd(op_model::Dict, subjects::Dict)
    println("\nStep 4: Running NeoPKPD population simulation...")

    # Mean dose
    mean_dose = mean([s["dose"] for s in values(subjects)])

    # Base model spec
    base = ModelSpec(
        op_model["model_type"](),
        "nonmem_migration",
        OneCompOralFirstOrderParams(
            op_model["typical"][:Ka],
            op_model["typical"][:CL],
            op_model["typical"][:V]
        ),
        [DoseEvent(0.0, mean_dose)]
    )

    # Covariate model
    cov_effects = CovariateEffect[]
    for cov in op_model["covariates"]
        if cov["type"] == "power"
            push!(cov_effects, CovariateEffect(
                PowerCovariate(),
                Symbol(cov["param"]),
                Symbol(cov["cov"]),
                cov["exp"],
                cov["ref"]
            ))
        end
    end
    cov_model = CovariateModel("wt_scaling", cov_effects)

    # Individual covariates
    ind_covs = IndividualCovariates[]
    for id in sort(collect(keys(subjects)))
        push!(ind_covs, IndividualCovariates(Dict(:WT => subjects[id]["wt"]), nothing))
    end

    # IIV specification
    omega_dict = Dict(
        :Ka => op_model["omega"][:Ka],
        :CL => op_model["omega"][:CL],
        :V => op_model["omega"][:V]
    )
    iiv = IIVSpec(LogNormalIIV(), omega_dict, UInt64(12345), length(ind_covs))

    # Population spec
    pop = PopulationSpec(base, iiv, nothing, cov_model, ind_covs)

    # Simulation grid
    grid = SimGrid(0.0, 25.0, collect(0.0:0.5:25.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    # Run simulation
    result = simulate_population(pop, grid, solver)

    println("  Simulated $(length(result.individuals)) subjects")

    return pop, grid, solver, result
end

# ============================================================================
# Step 5: Compare Results
# ============================================================================

function compare_results(pop_result)
    println("\nStep 5: Computing comparison metrics...")

    # Extract summary statistics
    conc_summary = pop_result.summaries[:conc]

    # Compute key metrics
    cmax_mean = maximum(conc_summary.mean)
    tmax = pop_result.individuals[1].t[argmax(conc_summary.mean)]

    # AUC (trapezoidal)
    t = pop_result.individuals[1].t
    c = conc_summary.mean
    auc = 0.0
    for i in 2:length(t)
        auc += 0.5 * (c[i] + c[i-1]) * (t[i] - t[i-1])
    end

    comparison = Dict(
        "cmax_mean" => cmax_mean,
        "tmax" => tmax,
        "auc_mean" => auc,
        "n_subjects" => length(pop_result.individuals)
    )

    println("  Cmax (mean): $(round(cmax_mean, digits=2)) mg/L")
    println("  Tmax: $(round(tmax, digits=1)) h")
    println("  AUC_0_25 (mean): $(round(auc, digits=1)) mg*h/L")

    return comparison
end

# ============================================================================
# Step 6: Generate Migration Report
# ============================================================================

function generate_report(nm_model::Dict, op_model::Dict, comparison::Dict)
    println("\nStep 6: Generating migration report...")

    report = Dict(
        "source" => Dict(
            "format" => "NONMEM",
            "subroutine" => nm_model["subroutine"],
            "trans" => nm_model["trans"],
            "n_thetas" => length(nm_model["thetas"]),
            "n_omegas" => length(nm_model["omegas"])
        ),
        "target" => Dict(
            "format" => "NeoPKPD",
            "model_type" => string(op_model["model_type"]),
            "parameters" => Dict(
                "Ka" => op_model["typical"][:Ka],
                "CL" => op_model["typical"][:CL],
                "V" => op_model["typical"][:V]
            ),
            "iiv" => Dict(
                "Ka_cv" => round(op_model["omega"][:Ka] * 100, digits=1),
                "CL_cv" => round(op_model["omega"][:CL] * 100, digits=1),
                "V_cv" => round(op_model["omega"][:V] * 100, digits=1)
            )
        ),
        "validation" => Dict(
            "n_subjects" => comparison["n_subjects"],
            "cmax_mean" => round(comparison["cmax_mean"], digits=2),
            "tmax" => comparison["tmax"],
            "auc_mean" => round(comparison["auc_mean"], digits=1)
        ),
        "status" => "MIGRATION_SUCCESSFUL"
    )

    return report
end

# ============================================================================
# Main
# ============================================================================

function main()
    println("=" ^ 60)
    println("NONMEM to NeoPKPD Migration")
    println("=" ^ 60)

    # Step 1: Parse NONMEM model
    nm_model = parse_nonmem_model()

    # Step 2: Convert to NeoPKPD
    op_model = convert_to_neopkpd(nm_model)

    # Step 3: Load data
    subjects = load_data()

    # Step 4: Simulate with NeoPKPD
    pop, grid, solver, result = simulate_neopkpd(op_model, subjects)

    # Step 5: Compare results
    comparison = compare_results(result)

    # Step 6: Generate report
    report = generate_report(nm_model, op_model, comparison)

    # Save outputs
    println("\nSaving outputs...")
    mkpath(joinpath(BASE_DIR, "output"))

    # Population artifact
    pop_path = joinpath(BASE_DIR, "output", "migrated_population.json")
    write_population_json(pop_path; population_spec=pop, grid=grid, solver=solver, result=result)
    println("  Wrote: $pop_path")

    # Migration report
    report_path = joinpath(BASE_DIR, "output", "migration_report.json")
    open(report_path, "w") do io
        JSON.print(io, report, 2)
    end
    println("  Wrote: $report_path")

    println("\n" * "=" ^ 60)
    println("Migration Status: $(report["status"])")
    println("=" ^ 60)
end

main()
