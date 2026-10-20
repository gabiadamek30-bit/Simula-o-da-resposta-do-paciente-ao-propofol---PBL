using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function weights_fixed()
    # Deterministic cohort: representative adult range
    return [50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
end

function make_covariate_model()
    return CovariateModel(
        "wt_on_cl_v",
        [
            CovariateEffect(PowerCovariate(), :CL, :WT, 0.75, 70.0),
            CovariateEffect(PowerCovariate(), :V, :WT, 1.00, 70.0),
        ],
    )
end

function run_one_dose(dose_mg::Float64)
    base = ModelSpec(
        OneCompIVBolus(),
        "fih_iv_" * string(Int(dose_mg)) * "mg",
        OneCompIVBolusParams(10.0, 50.0), # typical CL, V
        [DoseEvent(0.0, dose_mg)],
    )

    cm = make_covariate_model()

    covs = IndividualCovariates[]
    for wt in weights_fixed()
        push!(covs, IndividualCovariates(Dict(:WT => wt), nothing))
    end

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.25, :V => 0.20), UInt64(20260101), length(covs))

    pop = PopulationSpec(base, iiv, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    # Exposure metrics computed on mean concentration curve
    s = res.summaries[:conc]
    c_mean = s.mean
    t = res.individuals[1].t

    metrics = Dict(
        "dose_mg" => dose_mg,
        "cmax_mean" => cmax(t, c_mean),
        "auc0_24_mean" => auc_trapezoid(t, c_mean),
    )

    return pop, grid, solver, res, metrics
end

function main()
    base_dir = "docs/examples/use_cases/fih_dose_exploration"
    mkpath(joinpath(base_dir, "output"))

    dose_levels = [10.0, 30.0, 100.0, 300.0]

    all_metrics = Any[]

    for d in dose_levels
        pop, grid, solver, res, metrics = run_one_dose(d)
        push!(all_metrics, metrics)

        out_art = joinpath(base_dir, "output", "pop_iv_" * string(Int(d)) * "mg.json")
        write_population_json(out_art; population_spec = pop, grid = grid, solver = solver, result = res)
        println("Wrote artifact: " * out_art)
    end

    metrics_path = joinpath(base_dir, "output", "metrics.json")
    open(metrics_path, "w") do io
        JSON.print(io, all_metrics, 2)
    end
    println("Wrote metrics: " * metrics_path)
end

main()
