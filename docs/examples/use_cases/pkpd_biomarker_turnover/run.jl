using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

# PK/PD parameters
const CL_TYP = 10.0   # L/h
const V_TYP = 50.0    # L
const KIN = 10.0      # baseline production rate
const KOUT = 0.1      # elimination rate constant
const R0 = KIN / KOUT # baseline response = 100
const IMAX = 0.9      # maximal inhibition
const IC50 = 5.0      # concentration for half-maximal inhibition

# Population settings
const N_SUBJECTS = 20
const IIV_CL = 0.25
const IIV_V = 0.20
const IOV_CL = 0.15   # IOV on CL when enabled
const SEED = UInt64(20260103)

# Threshold for time_below: 80% of baseline
const THRESHOLD_80PCT = 0.8 * R0

function make_base_spec(name::String, doses::Vector{DoseEvent})
    return ModelSpec(
        OneCompIVBolus(),
        name,
        OneCompIVBolusParams(CL_TYP, V_TYP),
        doses,
    )
end

function make_pd_spec()
    # IndirectResponseTurnoverParams(Kin, Kout, R0, Imax, IC50)
    # PDSpec(kind, name, params, input_observation, output_observation)
    return PDSpec(
        IndirectResponseTurnover(),
        "biomarker_turnover",
        IndirectResponseTurnoverParams(KIN, KOUT, R0, IMAX, IC50),
        :conc,      # PD input is PK concentration
        :response,  # PD output is response
    )
end

function make_iiv(n::Int)
    return IIVSpec(LogNormalIIV(), Dict(:CL => IIV_CL, :V => IIV_V), SEED, n)
end

function make_iov()
    # IOVSpec uses same RandomEffectKind (LogNormalIIV) and OccasionDefinition
    return IOVSpec(LogNormalIIV(), Dict(:CL => IOV_CL), SEED + 1, OccasionDefinition(:dose_times))
end

function run_scenario(
    scenario::String,
    doses::Vector{DoseEvent},
    t_end::Float64,
    with_iov::Bool,
)
    base = make_base_spec(scenario, doses)
    pd = make_pd_spec()
    iiv = make_iiv(N_SUBJECTS)

    # IOV on CL across occasions (each dose time > t0 starts new occasion)
    iov = with_iov ? make_iov() : nothing

    # covariates must be a Vector{IndividualCovariates}, even if empty
    empty_covs = IndividualCovariates[]
    pop = PopulationSpec(base, iiv, iov, nothing, empty_covs)

    grid = SimGrid(0.0, t_end, collect(0.0:0.5:t_end))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver; pd_spec=pd)

    # Compute metrics on mean response curve
    t = res.individuals[1].t
    r_mean = res.summaries[:response].mean

    metrics = Dict(
        "scenario" => scenario,
        "with_iov" => with_iov,
        "emin_mean" => emin(t, r_mean),
        "time_below_80pct_mean" => time_below(t, r_mean, THRESHOLD_80PCT),
        "suppression_auc_mean" => auc_above_baseline(t, r_mean, R0),
    )

    return pop, grid, solver, pd, res, metrics
end

function main()
    base_dir = "docs/examples/use_cases/pkpd_biomarker_turnover"
    mkpath(joinpath(base_dir, "output"))

    # QD: 100 mg at 0, 24, 48
    qd_doses = [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0), DoseEvent(48.0, 100.0)]
    t_end_qd = 72.0

    # BID: 50 mg at 0, 12, 24, 36, 48, 60
    bid_doses = [
        DoseEvent(0.0, 50.0),
        DoseEvent(12.0, 50.0),
        DoseEvent(24.0, 50.0),
        DoseEvent(36.0, 50.0),
        DoseEvent(48.0, 50.0),
        DoseEvent(60.0, 50.0),
    ]
    t_end_bid = 72.0

    scenarios = [
        ("qd_100mg", qd_doses, t_end_qd),
        ("bid_50mg", bid_doses, t_end_bid),
    ]

    all_metrics = Any[]

    for (name, doses, t_end) in scenarios
        for with_iov in [false, true]
            suffix = with_iov ? "iov" : "no_iov"
            full_name = "$(name)_$(suffix)"

            println("Running: $(full_name)")
            pop, grid, solver, pd, res, metrics = run_scenario(name, doses, t_end, with_iov)
            push!(all_metrics, metrics)

            out_path = joinpath(base_dir, "output", "$(full_name).json")
            write_population_json(out_path; population_spec=pop, grid=grid, solver=solver, result=res, pd_spec=pd)
            println("  Wrote: $(out_path)")
        end
    end

    metrics_path = joinpath(base_dir, "output", "metrics.json")
    open(metrics_path, "w") do io
        JSON.print(io, all_metrics, 2)
    end
    println("Wrote metrics: $(metrics_path)")
end

main()
