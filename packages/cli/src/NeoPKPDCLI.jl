module NeoPKPDCLI

using ArgParse
using JSON
using NeoPKPD

# ============================================================================
# Utility Functions
# ============================================================================

function _die(msg::String)
    println(stderr, "Error: ", msg)
    exit(1)
end

function _info(msg::String)
    println(stderr, msg)
end

function _write_json(path::String, data::Any)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
end

# ============================================================================
# Version Command
# ============================================================================

const VERSION_HELP = """
Display version information for NeoPKPD.

Shows the NeoPKPD version along with semantic versions for event handling,
solver behavior, and artifact schema.

Example:
  neopkpd version
"""

function cmd_version()
    println("NeoPKPD " * NeoPKPD.NEOPKPD_VERSION)
    println("Event semantics: " * NeoPKPD.EVENT_SEMANTICS_VERSION)
    println("Solver semantics: " * NeoPKPD.SOLVER_SEMANTICS_VERSION)
    println("Artifact schema: " * NeoPKPD.ARTIFACT_SCHEMA_VERSION)
end

# ============================================================================
# Replay Command
# ============================================================================

const REPLAY_HELP = """
Replay a simulation from a saved artifact JSON file.

This command re-executes a simulation using the exact parameters stored in an
artifact file, validating reproducibility. Supports single, population, and
sensitivity artifacts.

Arguments:
  --artifact PATH    Path to the artifact JSON file (required)
  --out PATH         Optional output path for new artifact with replayed results

Examples:
  neopkpd replay --artifact validation/golden/pk_iv_bolus.json
  neopkpd replay --artifact my_sim.json --out replayed.json
"""

function cmd_replay(args)
    path = args["artifact"]
    out = get(args, "out", nothing)

    artifact = NeoPKPD.read_execution_json(path)
    atype = get(artifact, "artifact_type", "single")

    if atype == "population"
        res = NeoPKPD.replay_population_execution(artifact)
        _info("Replayed population simulation with $(length(res.individuals)) individuals")
        if out !== nothing
            parsed = NeoPKPD.deserialize_population_execution(artifact)
            NeoPKPD.write_population_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    if atype == "sensitivity_single"
        res = NeoPKPD.replay_sensitivity_execution(artifact)
        _info("Replayed single sensitivity analysis")
        if out !== nothing
            parsed = NeoPKPD.deserialize_sensitivity_execution(artifact)
            NeoPKPD.write_sensitivity_json(
                out;
                model_spec = parsed.model_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    if atype == "sensitivity_population"
        res = NeoPKPD.replay_population_sensitivity_execution(artifact)
        _info("Replayed population sensitivity analysis")
        if out !== nothing
            parsed = NeoPKPD.deserialize_population_sensitivity_execution(artifact)
            NeoPKPD.write_population_sensitivity_json(
                out;
                population_spec = parsed.population_spec,
                grid = parsed.grid,
                solver = parsed.solver,
                result = res,
            )
            _info("Written to: $out")
        end
        return
    end

    # default: single execution
    res = NeoPKPD.replay_execution(artifact)
    _info("Replayed single simulation")
    if out !== nothing
        parsed = NeoPKPD.deserialize_execution(artifact)
        NeoPKPD.write_execution_json(
            out;
            model_spec = parsed.model_spec,
            grid = parsed.grid,
            solver = parsed.solver,
            result = res,
            pd_spec = parsed.pd_spec,
        )
        _info("Written to: $out")
    end
end

# ============================================================================
# Validate Golden Command
# ============================================================================

const VALIDATE_GOLDEN_HELP = """
Run golden validation tests.

Validates all golden artifacts in the validation/golden directory to ensure
reproducibility across versions and platforms.

Example:
  neopkpd validate-golden
"""

function cmd_validate_golden()
    cli_src = @__DIR__
    repo_root = normpath(joinpath(cli_src, "..", "..", ".."))

    runner = joinpath(repo_root, "validation", "scripts", "run_golden_validation.jl")

    if !isfile(runner)
        _die("Golden validation runner not found: " * runner)
    end

    cmd = `julia $runner`
    run(Cmd(cmd; dir=repo_root))
end

# ============================================================================
# Simulate Command
# ============================================================================

const SIMULATE_HELP = """
Run a PK or PKPD simulation from a JSON specification file.

SUPPORTED PK MODELS:
  OneCompIVBolus          - One-compartment IV bolus
    params: {CL, V}
  OneCompOralFirstOrder   - One-compartment oral first-order absorption
    params: {Ka, CL, V}
  TwoCompIVBolus          - Two-compartment IV bolus
    params: {CL, V1, Q, V2}
  TwoCompOral             - Two-compartment oral first-order absorption
    params: {Ka, CL, V1, Q, V2}
  ThreeCompIVBolus        - Three-compartment IV bolus (mammillary)
    params: {CL, V1, Q2, V2, Q3, V3}
  TransitAbsorption       - Transit compartment chain model
    params: {N, Ktr, Ka, CL, V}
  MichaelisMentenElimination - Saturable (nonlinear) elimination
    params: {Vmax, Km, V}

SUPPORTED PD MODELS:
  DirectEmax              - Direct Emax (hyperbolic)
    params: {E0, Emax, EC50}
  SigmoidEmax             - Sigmoid Emax (Hill equation)
    params: {E0, Emax, EC50, gamma}
  IndirectResponseTurnover - Indirect response turnover
    params: {Kin, Kout, R0, Imax, IC50}
  BiophaseEquilibration   - Effect compartment model
    params: {ke0, E0, Emax, EC50}

Input JSON format:
{
  "model": {
    "kind": "<model_kind>",
    "name": "my_simulation",
    "params": {...},
    "doses": [{"time": 0.0, "amount": 100.0}]
  },
  "grid": {"t0": 0.0, "t1": 24.0, "saveat": [0.0, 1.0, 2.0, ...]},
  "solver": {"alg": "Tsit5", "reltol": 1e-10, "abstol": 1e-12, "maxiters": 10000000},
  "pd": {  # optional
    "kind": "<pd_kind>",
    "name": "pd_model",
    "params": {...},
    "input_observation": "conc",
    "output_observation": "effect"
  }
}

Arguments:
  --spec PATH       Path to simulation specification JSON (required)
  --out PATH        Output path for results JSON (required)
  --format FORMAT   Output format: "artifact" or "simple" (default: artifact)

Examples:
  neopkpd simulate --spec pk_spec.json --out result.json
  neopkpd simulate --spec pkpd_spec.json --out result.json --format simple

  # Two-compartment IV bolus example spec:
  # {"model": {"kind": "TwoCompIVBolus", "params": {"CL": 10, "V1": 50, "Q": 5, "V2": 100}, "doses": [{"time": 0, "amount": 500}]}, ...}

  # Sigmoid Emax PD example:
  # {"pd": {"kind": "SigmoidEmax", "params": {"E0": 0, "Emax": 100, "EC50": 5, "gamma": 2}, ...}}
"""

function _parse_model_spec(spec::Dict)
    kind_str = spec["kind"]
    name = get(spec, "name", "simulation")
    params_dict = spec["params"]
    doses_raw = get(spec, "doses", [])

    # Support both bolus and infusion doses - duration is optional (0.0 = bolus)
    doses = [NeoPKPD.DoseEvent(
        Float64(d["time"]),
        Float64(d["amount"]),
        Float64(get(d, "duration", 0.0))
    ) for d in doses_raw]

    if kind_str == "OneCompIVBolus"
        params = NeoPKPD.OneCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.OneCompIVBolus(), name, params, doses)

    elseif kind_str == "OneCompOralFirstOrder"
        params = NeoPKPD.OneCompOralFirstOrderParams(
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.OneCompOralFirstOrder(), name, params, doses)

    elseif kind_str == "TwoCompIVBolus"
        params = NeoPKPD.TwoCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q"]),
            Float64(params_dict["V2"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.TwoCompIVBolus(), name, params, doses)

    elseif kind_str == "TwoCompOral"
        params = NeoPKPD.TwoCompOralParams(
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q"]),
            Float64(params_dict["V2"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.TwoCompOral(), name, params, doses)

    elseif kind_str == "ThreeCompIVBolus"
        params = NeoPKPD.ThreeCompIVBolusParams(
            Float64(params_dict["CL"]),
            Float64(params_dict["V1"]),
            Float64(params_dict["Q2"]),
            Float64(params_dict["V2"]),
            Float64(params_dict["Q3"]),
            Float64(params_dict["V3"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.ThreeCompIVBolus(), name, params, doses)

    elseif kind_str == "TransitAbsorption"
        params = NeoPKPD.TransitAbsorptionParams(
            Int(params_dict["N"]),
            Float64(params_dict["Ktr"]),
            Float64(params_dict["Ka"]),
            Float64(params_dict["CL"]),
            Float64(params_dict["V"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.TransitAbsorption(), name, params, doses)

    elseif kind_str == "MichaelisMentenElimination"
        params = NeoPKPD.MichaelisMentenEliminationParams(
            Float64(params_dict["Vmax"]),
            Float64(params_dict["Km"]),
            Float64(params_dict["V"]),
        )
        return NeoPKPD.ModelSpec(NeoPKPD.MichaelisMentenElimination(), name, params, doses)

    else
        _die("Unsupported model kind: $kind_str. Supported: OneCompIVBolus, OneCompOralFirstOrder, TwoCompIVBolus, TwoCompOral, ThreeCompIVBolus, TransitAbsorption, MichaelisMentenElimination")
    end
end

function _parse_pd_spec(spec::Dict)
    kind_str = spec["kind"]
    name = get(spec, "name", "pd_model")
    params_dict = spec["params"]
    input_obs = Symbol(get(spec, "input_observation", "conc"))
    output_obs = Symbol(get(spec, "output_observation", "effect"))

    if kind_str == "DirectEmax"
        params = NeoPKPD.DirectEmaxParams(
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
        )
        return NeoPKPD.PDSpec(NeoPKPD.DirectEmax(), name, params, input_obs, output_obs)

    elseif kind_str == "SigmoidEmax"
        params = NeoPKPD.SigmoidEmaxParams(
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
            Float64(params_dict["gamma"]),
        )
        return NeoPKPD.PDSpec(NeoPKPD.SigmoidEmax(), name, params, input_obs, output_obs)

    elseif kind_str == "IndirectResponseTurnover"
        params = NeoPKPD.IndirectResponseTurnoverParams(
            Float64(params_dict["Kin"]),
            Float64(params_dict["Kout"]),
            Float64(params_dict["R0"]),
            Float64(params_dict["Imax"]),
            Float64(params_dict["IC50"]),
        )
        return NeoPKPD.PDSpec(NeoPKPD.IndirectResponseTurnover(), name, params, input_obs, output_obs)

    elseif kind_str == "BiophaseEquilibration"
        params = NeoPKPD.BiophaseEquilibrationParams(
            Float64(params_dict["ke0"]),
            Float64(params_dict["E0"]),
            Float64(params_dict["Emax"]),
            Float64(params_dict["EC50"]),
        )
        return NeoPKPD.PDSpec(NeoPKPD.BiophaseEquilibration(), name, params, input_obs, output_obs)

    else
        _die("Unsupported PD kind: $kind_str. Supported: DirectEmax, SigmoidEmax, IndirectResponseTurnover, BiophaseEquilibration")
    end
end

function _parse_grid(spec::Dict)
    return NeoPKPD.SimGrid(
        Float64(spec["t0"]),
        Float64(spec["t1"]),
        [Float64(x) for x in spec["saveat"]],
    )
end

function _parse_solver(spec::Dict)
    return NeoPKPD.SolverSpec(
        Symbol(get(spec, "alg", "Tsit5")),
        Float64(get(spec, "reltol", 1e-10)),
        Float64(get(spec, "abstol", 1e-12)),
        Int(get(spec, "maxiters", 10^7)),
    )
end

function _simresult_to_dict(res::NeoPKPD.SimResult)
    return Dict(
        "t" => res.t,
        "states" => Dict(string(k) => v for (k, v) in res.states),
        "observations" => Dict(string(k) => v for (k, v) in res.observations),
        "metadata" => res.metadata,
    )
end

function cmd_simulate(args)
    spec_path = args["spec"]
    out_path = args["out"]
    format = get(args, "format", "artifact")

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : NeoPKPD.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd_spec = haskey(spec, "pd") ? _parse_pd_spec(spec["pd"]) : nothing

    if pd_spec !== nothing
        pd_kind = pd_spec.kind
        # Determine if this is a direct (pk_then_pd) or coupled PD model
        if pd_kind isa NeoPKPD.DirectEmax || pd_kind isa NeoPKPD.SigmoidEmax || pd_kind isa NeoPKPD.BiophaseEquilibration
            # Direct PD models - use simulate_pkpd (pk_then_pd mode)
            result = NeoPKPD.simulate_pkpd(model_spec, pd_spec, grid, solver)
            _info("Simulated PK then PD: $(length(result.t)) time points")
        elseif pd_kind isa NeoPKPD.IndirectResponseTurnover
            # Coupled ODE PD models - use simulate_pkpd_coupled
            result = NeoPKPD.simulate_pkpd_coupled(model_spec, pd_spec, grid, solver)
            _info("Simulated coupled PKPD: $(length(result.t)) time points")
        else
            _die("Unsupported PD model type for simulation")
        end
    else
        # PK only
        result = NeoPKPD.simulate(model_spec, grid, solver)
        _info("Simulated PK: $(length(result.t)) time points")
    end

    if format == "simple"
        _write_json(out_path, _simresult_to_dict(result))
    else
        NeoPKPD.write_execution_json(
            out_path;
            model_spec = model_spec,
            grid = grid,
            solver = solver,
            result = result,
            pd_spec = pd_spec,
        )
    end

    _info("Written to: $out_path")
end

# ============================================================================
# Population Command
# ============================================================================

const POPULATION_HELP = """
Run a population PK/PD simulation with inter-individual variability.

Input JSON format:
{
  "model": { ... },  # Same as simulate command
  "grid": { ... },
  "solver": { ... },
  "iiv": {
    "kind": "LogNormalIIV",
    "omegas": {"CL": 0.3, "V": 0.2},
    "seed": 12345,
    "n": 100
  },
  "iov": {  # optional
    "kind": "LogNormalIIV",
    "pis": {"CL": 0.1},
    "seed": 54321,
    "occasion_def": "dose_times"
  },
  "covariate_model": {  # optional
    "name": "wt_model",
    "effects": [
      {"kind": "PowerCovariate", "param": "CL", "covariate": "WT", "beta": 0.75, "ref": 70.0}
    ]
  },
  "covariates": [  # optional, one per individual
    {"values": {"WT": 70.0}},
    ...
  ]
}

Arguments:
  --spec PATH       Path to population specification JSON (required)
  --out PATH        Output path for results JSON (required)
  --format FORMAT   Output format: "artifact" or "simple" (default: artifact)

Examples:
  neopkpd population --spec pop_spec.json --out pop_result.json
"""

function _parse_iiv_spec(spec::Dict)
    omegas = Dict(Symbol(k) => Float64(v) for (k, v) in spec["omegas"])
    return NeoPKPD.IIVSpec(
        NeoPKPD.LogNormalIIV(),
        omegas,
        UInt64(spec["seed"]),
        Int(spec["n"]),
    )
end

function _parse_iov_spec(spec::Dict)
    pis = Dict(Symbol(k) => Float64(v) for (k, v) in spec["pis"])
    occasion_mode = Symbol(get(spec, "occasion_def", "dose_times"))
    return NeoPKPD.IOVSpec(
        NeoPKPD.LogNormalIIV(),
        pis,
        UInt64(spec["seed"]),
        NeoPKPD.OccasionDefinition(occasion_mode),
    )
end

function _parse_covariate_effect(spec::Dict)
    kind_str = spec["kind"]
    if kind_str == "LinearCovariate"
        kind = NeoPKPD.LinearCovariate()
    elseif kind_str == "PowerCovariate"
        kind = NeoPKPD.PowerCovariate()
    elseif kind_str == "ExpCovariate"
        kind = NeoPKPD.ExpCovariate()
    else
        _die("Unsupported covariate effect kind: $kind_str")
    end
    return NeoPKPD.CovariateEffect(
        kind,
        Symbol(spec["param"]),
        Symbol(spec["covariate"]),
        Float64(spec["beta"]),
        Float64(spec["ref"]),
    )
end

function _parse_covariate_model(spec::Dict)
    effects = [_parse_covariate_effect(e) for e in spec["effects"]]
    return NeoPKPD.CovariateModel(spec["name"], effects)
end

function _parse_individual_covariates(specs::Vector)
    return [
        NeoPKPD.IndividualCovariates(
            Dict(Symbol(k) => Float64(v) for (k, v) in get(s, "values", Dict())),
            nothing,  # time_varying not supported via CLI yet
        )
        for s in specs
    ]
end

function _popresult_to_dict(res)
    individuals = [_simresult_to_dict(r) for r in res.individuals]
    params = [Dict(string(k) => v for (k, v) in d) for d in res.params]

    summaries = Dict()
    for (k, s) in res.summaries
        summaries[string(k)] = Dict(
            "observation" => string(s.observation),
            "probs" => s.probs,
            "mean" => s.mean,
            "median" => s.median,
            "quantiles" => Dict(string(p) => v for (p, v) in s.quantiles),
        )
    end

    return Dict(
        "individuals" => individuals,
        "params" => params,
        "summaries" => summaries,
        "metadata" => res.metadata,
    )
end

function cmd_population(args)
    spec_path = args["spec"]
    out_path = args["out"]
    format = get(args, "format", "artifact")

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : NeoPKPD.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    iiv = haskey(spec, "iiv") ? _parse_iiv_spec(spec["iiv"]) : nothing
    iov = haskey(spec, "iov") ? _parse_iov_spec(spec["iov"]) : nothing
    cov_model = haskey(spec, "covariate_model") ? _parse_covariate_model(spec["covariate_model"]) : nothing
    covariates = haskey(spec, "covariates") ? _parse_individual_covariates(spec["covariates"]) : NeoPKPD.IndividualCovariates[]

    pop_spec = NeoPKPD.PopulationSpec(model_spec, iiv, iov, cov_model, covariates)

    result = NeoPKPD.simulate_population(pop_spec, grid, solver)
    n = length(result.individuals)
    _info("Simulated population: $n individuals")

    if format == "simple"
        _write_json(out_path, _popresult_to_dict(result))
    else
        NeoPKPD.write_population_json(
            out_path;
            population_spec = pop_spec,
            grid = grid,
            solver = solver,
            result = result,
        )
    end

    _info("Written to: $out_path")
end

# ============================================================================
# Sensitivity Command
# ============================================================================

const SENSITIVITY_HELP = """
Run sensitivity analysis by perturbing parameters.

Input JSON format:
{
  "model": { ... },  # Same as simulate command
  "grid": { ... },
  "solver": { ... },
  "perturbation": {
    "name": "cl_sensitivity",
    "param": "CL",
    "delta": 0.01
  },
  "observation": "conc"
}

Arguments:
  --spec PATH       Path to sensitivity specification JSON (required)
  --out PATH        Output path for results JSON (required)

Examples:
  neopkpd sensitivity --spec sens_spec.json --out sens_result.json
"""

function cmd_sensitivity(args)
    spec_path = args["spec"]
    out_path = args["out"]

    if !isfile(spec_path)
        _die("Specification file not found: $spec_path")
    end

    spec = JSON.parsefile(spec_path; dicttype=Dict{String, Any})

    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : NeoPKPD.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pert_spec = spec["perturbation"]
    plan = NeoPKPD.PerturbationPlan(
        pert_spec["name"],
        Symbol(pert_spec["param"]),
        Float64(pert_spec["delta"]),
    )

    observation = Symbol(get(spec, "observation", "conc"))

    result = NeoPKPD.run_sensitivity(model_spec, grid, solver, plan, observation)

    _info("Sensitivity analysis complete")
    _info("  Max absolute delta: $(result.metrics.max_abs_delta)")
    _info("  Max relative delta: $(result.metrics.max_rel_delta)")
    _info("  L2 norm delta: $(result.metrics.l2_norm_delta)")

    NeoPKPD.write_sensitivity_json(
        out_path;
        model_spec = model_spec,
        grid = grid,
        solver = solver,
        result = result,
    )

    _info("Written to: $out_path")
end

# ============================================================================
# Metrics Command
# ============================================================================

const METRICS_HELP = """
Compute PK/PD metrics from a simulation result or artifact.

Available metrics:
  - cmax: Maximum concentration
  - auc: Area under the curve (trapezoidal)
  - emin: Minimum response value
  - time_below: Time spent below a threshold
  - auc_above_baseline: AUC above a baseline value

Arguments:
  --artifact PATH      Path to artifact or result JSON (required)
  --observation NAME   Observation to analyze (default: conc)
  --metric NAME        Metric to compute: cmax, auc, emin, time_below, auc_above_baseline
  --threshold VALUE    Threshold for time_below or auc_above_baseline

Examples:
  neopkpd metrics --artifact result.json --metric cmax
  neopkpd metrics --artifact result.json --observation effect --metric emin
  neopkpd metrics --artifact result.json --metric time_below --threshold 10.0
"""

function cmd_metrics(args)
    artifact_path = args["artifact"]
    observation = Symbol(get(args, "observation", "conc"))
    metric_name = get(args, "metric", nothing)
    threshold = haskey(args, "threshold") && args["threshold"] !== nothing ? parse(Float64, args["threshold"]) : nothing

    if !isfile(artifact_path)
        _die("Artifact file not found: $artifact_path")
    end

    if metric_name === nothing
        _die("Metric name required. Use --metric cmax|auc|emin|time_below|auc_above_baseline")
    end

    artifact = JSON.parsefile(artifact_path)

    # Extract time and observation data
    t = nothing
    y = nothing

    if haskey(artifact, "result") && haskey(artifact["result"], "t")
        # Artifact format
        t = Float64.(artifact["result"]["t"])
        obs_data = artifact["result"]["observations"]
        obs_key = string(observation)
        if !haskey(obs_data, obs_key)
            _die("Observation '$obs_key' not found. Available: $(join(keys(obs_data), ", "))")
        end
        y = Float64.(obs_data[obs_key])
    elseif haskey(artifact, "t")
        # Simple format
        t = Float64.(artifact["t"])
        obs_data = artifact["observations"]
        obs_key = string(observation)
        if !haskey(obs_data, obs_key)
            _die("Observation '$obs_key' not found. Available: $(join(keys(obs_data), ", "))")
        end
        y = Float64.(obs_data[obs_key])
    else
        _die("Unrecognized file format")
    end

    result = nothing

    if metric_name == "cmax"
        result = NeoPKPD.cmax(t, y)
        println("Cmax: $result")
    elseif metric_name == "auc"
        result = NeoPKPD.auc_trapezoid(t, y)
        println("AUC: $result")
    elseif metric_name == "emin"
        result = NeoPKPD.emin(t, y)
        println("Emin: $result")
    elseif metric_name == "time_below"
        if threshold === nothing
            _die("--threshold required for time_below metric")
        end
        result = NeoPKPD.time_below(t, y, threshold)
        println("Time below $threshold: $result")
    elseif metric_name == "auc_above_baseline"
        if threshold === nothing
            _die("--threshold required for auc_above_baseline metric")
        end
        result = NeoPKPD.auc_above_baseline(t, y, threshold)
        println("AUC above baseline $threshold: $result")
    else
        _die("Unknown metric: $metric_name. Supported: cmax, auc, emin, time_below, auc_above_baseline")
    end
end

# ============================================================================
# Import Command (NONMEM / Monolix)
# ============================================================================

const IMPORT_HELP = """
Import a model from NONMEM control file or Monolix project.

This command parses external model files and converts them to NeoPKPD format.

Arguments:
  --file PATH      Path to .ctl (NONMEM) or .mlxtran (Monolix) file (required)
  --format FORMAT  Input format: nonmem or monolix (auto-detected if not specified)
  --out PATH       Output path for NeoPKPD JSON spec (required)

Examples:
  neopkpd import --file model.ctl --out neopkpd_spec.json
  neopkpd import --file project.mlxtran --format monolix --out spec.json
"""

function cmd_import(args)
    filepath = args["file"]
    out = args["out"]
    format = get(args, "format", nothing)

    if !isfile(filepath)
        _die("File not found: $filepath")
    end

    # Auto-detect format
    if format === nothing
        if endswith(filepath, ".ctl") || endswith(filepath, ".mod")
            format = "nonmem"
        elseif endswith(filepath, ".mlxtran")
            format = "monolix"
        else
            _die("Cannot auto-detect format. Please specify --format nonmem or --format monolix")
        end
    end

    if format == "nonmem"
        result = NeoPKPD.parse_nonmem_control(filepath)
        model_spec, pop_spec, metadata = NeoPKPD.convert_nonmem_to_neopkpd(result)
        _info("Parsed NONMEM control file: $(result.problem)")
    elseif format == "monolix"
        result = NeoPKPD.parse_monolix_project(filepath)
        model_spec, pop_spec, metadata = NeoPKPD.convert_monolix_to_neopkpd(result)
        _info("Parsed Monolix project")
    else
        _die("Unknown format: $format. Use 'nonmem' or 'monolix'")
    end

    # Write output
    output = Dict(
        "model" => Dict(
            "kind" => string(typeof(model_spec.kind)),
            "params" => model_spec.params,
            "doses" => model_spec.doses
        ),
        "metadata" => metadata
    )

    if pop_spec !== nothing
        output["population"] = Dict(
            "iiv" => pop_spec.iiv !== nothing ? Dict(
                "kind" => string(typeof(pop_spec.iiv.kind)),
                "omegas" => pop_spec.iiv.omegas,
                "n" => pop_spec.iiv.n
            ) : nothing
        )
    end

    _write_json(out, output)
    _info("Written NeoPKPD specification to: $out")
end

# ============================================================================
# NCA Command (Non-Compartmental Analysis)
# ============================================================================

const NCA_HELP = """
Run Non-Compartmental Analysis on observed or simulated data.

Arguments:
  --data PATH       Path to data file (CSV with time, conc columns) (required)
  --dose FLOAT      Dose amount for CL/F and Vz/F calculation (required)
  --out PATH        Output path for NCA results (required)
  --method METHOD   AUC method: linear, log_linear, linear_up_log_down (default: linear_up_log_down)
  --tau FLOAT       Dosing interval for steady-state metrics (optional)

Examples:
  neopkpd nca --data pk_data.csv --dose 100 --out nca_results.json
  neopkpd nca --data ss_data.csv --dose 100 --tau 24 --out nca_ss.json
"""

function cmd_nca(args)
    data_path = args["data"]
    dose = parse(Float64, args["dose"])
    out = args["out"]
    method_str = get(args, "method", "linear_up_log_down")
    tau_str = get(args, "tau", nothing)

    if !isfile(data_path)
        _die("Data file not found: $data_path")
    end

    # Read CSV data
    lines = readlines(data_path)
    if length(lines) < 2
        _die("Data file must have header and at least one data row")
    end

    header = split(lines[1], ",")
    time_idx = findfirst(x -> lowercase(strip(x)) in ["time", "t"], header)
    conc_idx = findfirst(x -> lowercase(strip(x)) in ["conc", "concentration", "dv", "y"], header)

    if time_idx === nothing || conc_idx === nothing
        _die("Data file must have 'time' and 'conc' columns")
    end

    times = Float64[]
    concs = Float64[]
    for line in lines[2:end]
        parts = split(line, ",")
        if length(parts) >= max(time_idx, conc_idx)
            t = tryparse(Float64, strip(parts[time_idx]))
            c = tryparse(Float64, strip(parts[conc_idx]))
            if t !== nothing && c !== nothing
                push!(times, t)
                push!(concs, c)
            end
        end
    end

    if length(times) < 3
        _die("Need at least 3 data points for NCA")
    end

    # Set up NCA method
    method = if method_str == "linear"
        NeoPKPD.LinearMethod()
    elseif method_str == "log_linear"
        NeoPKPD.LogLinearMethod()
    else
        NeoPKPD.LinearUpLogDownMethod()
    end

    # Run NCA
    tau = tau_str !== nothing ? parse(Float64, tau_str) : nothing

    nca_result = NeoPKPD.run_nca(times, concs, dose; method=method, tau=tau)

    # Build output
    output = Dict(
        "cmax" => nca_result.cmax,
        "tmax" => nca_result.tmax,
        "clast" => nca_result.clast,
        "tlast" => nca_result.tlast,
        "auc_0_t" => nca_result.auc_0_t,
        "auc_0_inf" => nca_result.auc_0_inf,
        "aumc_0_t" => nca_result.aumc_0_t,
        "aumc_0_inf" => nca_result.aumc_0_inf,
        "lambda_z" => nca_result.lambda_z,
        "half_life" => nca_result.half_life,
        "mrt" => nca_result.mrt,
        "cl_f" => nca_result.cl_f,
        "vz_f" => nca_result.vz_f,
        "vss" => nca_result.vss,
        "method" => method_str,
        "dose" => dose,
        "n_points" => length(times)
    )

    if tau !== nothing
        output["tau"] = tau
        output["auc_0_tau"] = get(nca_result, :auc_0_tau, nothing)
        output["cavg"] = get(nca_result, :cavg, nothing)
    end

    _write_json(out, output)
    _info("NCA results written to: $out")
    println("Cmax: $(nca_result.cmax)")
    println("AUC(0-inf): $(nca_result.auc_0_inf)")
    println("Half-life: $(nca_result.half_life)")
    println("CL/F: $(nca_result.cl_f)")
end

# ============================================================================
# Helper: Parse Population Specification
# ============================================================================

function _parse_population_spec(spec::Dict)
    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : NeoPKPD.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    iiv = haskey(spec, "iiv") ? _parse_iiv_spec(spec["iiv"]) : nothing
    iov = haskey(spec, "iov") ? _parse_iov_spec(spec["iov"]) : nothing
    cov_model = haskey(spec, "covariate_model") ? _parse_covariate_model(spec["covariate_model"]) : nothing
    covariates = haskey(spec, "covariates") ? _parse_individual_covariates(spec["covariates"]) : NeoPKPD.IndividualCovariates[]

    pop_spec = NeoPKPD.PopulationSpec(model_spec, iiv, iov, cov_model, covariates)

    return pop_spec, grid, solver
end

function _parse_estimation_model_spec(spec::Dict)
    model_spec = _parse_model_spec(spec["model"])
    grid = _parse_grid(spec["grid"])
    solver = haskey(spec, "solver") ? _parse_solver(spec["solver"]) : NeoPKPD.SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)
    return model_spec, grid, solver
end

# ============================================================================
# VPC Command (Visual Predictive Check)
# ============================================================================

const VPC_HELP = """
Run Visual Predictive Check analysis.

Arguments:
  --observed PATH     Path to observed data CSV (required)
  --pop-spec PATH     Path to population specification JSON (required)
  --out PATH          Output path for VPC results (required)
  --n-sim INT         Number of simulations (default: 200)
  --n-bins INT        Number of time bins (default: 8)
  --pi LEVELS         Prediction interval levels, comma-separated (default: 0.05,0.50,0.95)
  --pcvpc             Use prediction-corrected VPC

Examples:
  neopkpd vpc --observed data.csv --pop-spec pop.json --out vpc.json
  neopkpd vpc --observed data.csv --pop-spec pop.json --out vpc.json --pcvpc --n-sim 500
"""

function cmd_vpc(args)
    observed_path = args["observed"]
    pop_spec_path = args["pop-spec"]
    out = args["out"]
    n_sim = parse(Int, get(args, "n-sim", "200"))
    n_bins = parse(Int, get(args, "n-bins", "8"))
    pi_str = get(args, "pi", "0.05,0.50,0.95")
    pcvpc = get(args, "pcvpc", false)

    if !isfile(observed_path)
        _die("Observed data file not found: $observed_path")
    end
    if !isfile(pop_spec_path)
        _die("Population spec file not found: $pop_spec_path")
    end

    # Parse PI levels
    pi_levels = [parse(Float64, strip(x)) for x in split(pi_str, ",")]

    # Read observed data
    observed = _read_observed_data(observed_path)

    # Read population spec
    pop_spec_json = JSON.parsefile(pop_spec_path)
    pop_spec, grid, solver = _parse_population_spec(pop_spec_json)

    # VPC config
    config = NeoPKPD.VPCConfig(
        pi_levels = pi_levels,
        binning = NeoPKPD.QuantileBinning(n_bins),
        n_simulations = n_sim,
        n_bootstrap = 500,
        prediction_corrected = pcvpc,
        seed = UInt64(12345)
    )

    # Run VPC
    vpc_result = if pcvpc
        NeoPKPD.compute_pcvpc(observed, pop_spec, grid, solver; config=config)
    else
        NeoPKPD.compute_vpc(observed, pop_spec, grid, solver; config=config)
    end

    # Build output
    output = Dict(
        "n_subjects" => vpc_result.n_subjects_observed,
        "n_observations" => vpc_result.n_observations_observed,
        "n_simulations" => vpc_result.n_simulations,
        "prediction_corrected" => pcvpc,
        "pi_levels" => pi_levels,
        "bins" => [
            Dict(
                "id" => bin.bin_id,
                "time_min" => bin.time_min,
                "time_max" => bin.time_max,
                "time_midpoint" => bin.time_midpoint,
                "n_observed" => bin.n_observed,
                "percentiles" => [
                    Dict(
                        "level" => p.level,
                        "observed" => p.observed,
                        "simulated_median" => p.simulated_median,
                        "simulated_lower" => p.simulated_lower,
                        "simulated_upper" => p.simulated_upper
                    ) for p in bin.percentiles
                ]
            ) for bin in vpc_result.bins
        ]
    )

    _write_json(out, output)
    _info("VPC results written to: $out")
    _info("Bins: $(length(vpc_result.bins)), Subjects: $(vpc_result.n_subjects_observed)")
end

function _read_observed_data(path::String)
    lines = readlines(path)
    header = split(lines[1], ",")

    id_idx = findfirst(x -> lowercase(strip(x)) in ["id", "subject", "usubjid"], header)
    time_idx = findfirst(x -> lowercase(strip(x)) in ["time", "t"], header)
    dv_idx = findfirst(x -> lowercase(strip(x)) in ["dv", "conc", "y", "observation"], header)

    if time_idx === nothing || dv_idx === nothing
        error("Data must have 'time' and 'dv' (or 'conc') columns")
    end

    # Group by subject
    subjects = Dict{String, Tuple{Vector{Float64}, Vector{Float64}}}()
    default_doses = [NeoPKPD.DoseEvent(0.0, 100.0)]

    for line in lines[2:end]
        parts = split(line, ",")
        if length(parts) >= max(time_idx, dv_idx)
            subj_id = id_idx !== nothing ? strip(parts[id_idx]) : "SUBJ001"
            t = tryparse(Float64, strip(parts[time_idx]))
            dv = tryparse(Float64, strip(parts[dv_idx]))

            if t !== nothing && dv !== nothing
                if !haskey(subjects, subj_id)
                    subjects[subj_id] = (Float64[], Float64[])
                end
                push!(subjects[subj_id][1], t)
                push!(subjects[subj_id][2], dv)
            end
        end
    end

    # Convert to ObservedData
    subject_data = [
        NeoPKPD.SubjectData(id, times, obs, default_doses)
        for (id, (times, obs)) in subjects
    ]

    return NeoPKPD.ObservedData(subject_data)
end

# ============================================================================
# Estimate Command (Parameter Estimation)
# ============================================================================

const ESTIMATE_HELP = """
Run parameter estimation using NLME methods.

Arguments:
  --data PATH         Path to observed data CSV (required)
  --model-spec PATH   Path to model specification JSON (required)
  --out PATH          Output path for estimation results (required)
  --method METHOD     Estimation method: laplacian, foce, saem (default: foce)
  --max-iter INT      Maximum iterations (default: 500)
  --compute-se        Compute standard errors (default: true)

Examples:
  neopkpd estimate --data pk_data.csv --model-spec model.json --out est.json
  neopkpd estimate --data data.csv --model-spec model.json --out est.json --method saem
"""

function cmd_estimate(args)
    data_path = args["data"]
    model_spec_path = args["model-spec"]
    out = args["out"]
    method_str = get(args, "method", "foce")
    max_iter = parse(Int, get(args, "max-iter", "500"))
    compute_se = get(args, "compute-se", "true") == "true"

    if !isfile(data_path)
        _die("Data file not found: $data_path")
    end
    if !isfile(model_spec_path)
        _die("Model spec file not found: $model_spec_path")
    end

    # Read observed data
    observed = _read_observed_data(data_path)

    # Read model specification
    spec_json = JSON.parsefile(model_spec_path)
    model_spec, grid, solver = _parse_estimation_model_spec(spec_json)

    # Get initial values from spec
    theta_init = get(spec_json, "theta_init", nothing)
    omega_init = get(spec_json, "omega_init", nothing)

    if theta_init === nothing
        _die("Model spec must include 'theta_init' for initial parameter values")
    end
    if omega_init === nothing
        _die("Model spec must include 'omega_init' for initial omega matrix")
    end

    theta_init = Float64.(theta_init)
    omega_init = Float64.(hcat(omega_init...))

    # Set up estimation method
    method = if method_str == "laplacian"
        NeoPKPD.LaplacianMethod()
    elseif method_str == "saem"
        NeoPKPD.SAEMMethod(n_burn=100, n_iter=100)
    else
        NeoPKPD.FOCEIMethod()
    end

    # Set up sigma
    sigma_init = NeoPKPD.ResidualErrorSpec(
        NeoPKPD.ProportionalError(),
        NeoPKPD.ProportionalErrorParams(0.1),
        :conc,
        UInt64(1)
    )

    # Estimation config
    config = NeoPKPD.EstimationConfig(
        method;
        theta_init = theta_init,
        omega_init = omega_init,
        sigma_init = sigma_init,
        max_iter = max_iter,
        compute_se = compute_se,
        verbose = true
    )

    _info("Running estimation with $(typeof(method))...")
    result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

    # Build output
    output = Dict(
        "method" => method_str,
        "convergence" => result.convergence,
        "n_iterations" => result.n_iterations,
        "ofv" => result.ofv,
        "aic" => result.aic,
        "bic" => result.bic,
        "theta" => result.theta,
        "theta_se" => result.theta_se,
        "theta_rse" => result.theta_rse,
        "omega" => [collect(row) for row in eachrow(result.omega)],
        "omega_corr" => [collect(row) for row in eachrow(result.omega_corr)],
        "n_subjects" => length(result.individuals),
        "runtime_seconds" => result.runtime_seconds
    )

    if result.theta_ci_lower !== nothing
        output["theta_ci_lower"] = result.theta_ci_lower
        output["theta_ci_upper"] = result.theta_ci_upper
    end

    _write_json(out, output)
    _info("Estimation results written to: $out")
    println("Convergence: $(result.convergence)")
    println("OFV: $(round(result.ofv, digits=3))")
    println("Theta: $(round.(result.theta, digits=4))")
    if result.theta_rse !== nothing
        println("RSE%: $(round.(result.theta_rse, digits=1))")
    end
end

# ============================================================================
# Trial Simulation Command
# ============================================================================

const TRIAL_HELP = """
Run clinical trial simulation.

Arguments:
  --spec PATH         Path to trial specification JSON (required)
  --out PATH          Output path for trial results (required)
  --replicates N      Number of trial replicates for power analysis (default: 1)

Supported study designs:
  - ParallelDesign: Parallel group studies
  - CrossoverDesign: Crossover studies
  - DoseEscalationDesign: Phase I dose escalation (3+3, mTPI, CRM)
  - BioequivalenceDesign: BE studies

Trial spec structure:
  {
    "name": "Phase 2 Study",
    "design": {"type": "parallel", "n_arms": 2},
    "arms": [
      {
        "name": "Placebo",
        "dose": 0.0,
        "n_subjects": 50
      },
      {
        "name": "Active",
        "dose": 100.0,
        "n_subjects": 50
      }
    ],
    "duration_days": 28,
    "pk_sampling_times": [0, 1, 2, 4, 8, 12, 24],
    "endpoints": [{"name": "pk_auc", "type": "pk", "metric": "auc_0_inf"}]
  }

Examples:
  neopkpd trial --spec trial_spec.json --out trial_result.json
  neopkpd trial --spec trial_spec.json --out trial_result.json --replicates 100
"""

function cmd_trial(args)
    spec_path = args["spec"]
    out = args["out"]
    n_replicates = parse(Int, get(args, "replicates", "1"))

    if !isfile(spec_path)
        _die("Spec file not found: $spec_path")
    end

    spec_dict = JSON.parsefile(spec_path)
    _info("Running trial simulation: $(get(spec_dict, "name", "Unnamed Trial"))")

    # Parse trial specification
    design = _parse_trial_design(spec_dict["design"])
    arms = [_parse_treatment_arm(a) for a in spec_dict["arms"]]
    endpoints = haskey(spec_dict, "endpoints") ?
        [_parse_endpoint(e) for e in spec_dict["endpoints"]] :
        NeoPKPD.EndpointSpec[]

    duration_days = Float64(get(spec_dict, "duration_days", 28.0))
    pk_sampling = haskey(spec_dict, "pk_sampling_times") ?
        [Float64(t) for t in spec_dict["pk_sampling_times"]] :
        [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]

    trial_spec = NeoPKPD.TrialSpec(
        get(spec_dict, "name", "CLI Trial"),
        design,
        arms;
        duration_days = duration_days,
        pk_sampling_times = pk_sampling,
        endpoints = endpoints,
        n_replicates = n_replicates,
        seed = UInt64(get(spec_dict, "seed", 12345))
    )

    # Run trial simulation
    if n_replicates > 1
        _info("Running power simulation with $n_replicates replicates...")
        result = NeoPKPD.run_power_simulation(trial_spec)
    else
        result = NeoPKPD.simulate_trial(trial_spec)
    end

    # Build output
    output = Dict(
        "trial_name" => result.trial_name,
        "design_type" => result.design_type,
        "n_replicates" => result.n_replicates,
        "seed" => result.seed,
        "arms" => Dict(
            name => Dict(
                "n_enrolled" => arm.n_enrolled,
                "n_completed" => arm.n_completed,
                "n_dropout" => arm.n_dropout,
                "summary_stats" => arm.summary_stats
            ) for (name, arm) in result.arms
        ),
        "endpoint_analyses" => result.endpoint_analyses,
        "power_estimates" => result.power_estimates
    )

    if result.be_results !== nothing
        output["be_results"] = result.be_results
    end

    _write_json(out, output)
    _info("Trial results written to: $out")

    # Print summary
    for (arm_name, arm_result) in result.arms
        println("  $arm_name: $(arm_result.n_completed)/$(arm_result.n_enrolled) completed")
    end

    if !isempty(result.power_estimates)
        println("Power estimates:")
        for (endpoint, power) in result.power_estimates
            println("  $endpoint: $(round(power * 100, digits=1))%")
        end
    end
end

function _parse_trial_design(design_dict)
    design_type = lowercase(get(design_dict, "type", "parallel"))

    if design_type == "parallel"
        n_arms = get(design_dict, "n_arms", 2)
        return NeoPKPD.ParallelDesign(n_arms)
    elseif design_type == "crossover"
        n_periods = get(design_dict, "n_periods", 2)
        n_sequences = get(design_dict, "n_sequences", 2)
        washout = Float64(get(design_dict, "washout_duration", 7.0))
        return NeoPKPD.CrossoverDesign(n_periods, n_sequences; washout_duration=washout)
    elseif design_type == "dose_escalation"
        dose_levels = [Float64(d) for d in design_dict["dose_levels"]]
        return NeoPKPD.DoseEscalationDesign(dose_levels)
    elseif design_type == "bioequivalence"
        return NeoPKPD.BioequivalenceDesign()
    else
        error("Unknown design type: $design_type")
    end
end

function _parse_treatment_arm(arm_dict)
    name = arm_dict["name"]
    dose = Float64(get(arm_dict, "dose", 100.0))
    n_subjects = get(arm_dict, "n_subjects", 50)

    # Create a simple model spec for the arm
    params = NeoPKPD.OneCompIVBolusParams(10.0, 100.0)  # Default CL, V
    doses = [NeoPKPD.DoseEvent(0.0, dose)]
    model_spec = NeoPKPD.ModelSpec(NeoPKPD.OneCompIVBolus(), name, params, doses)

    regimen = NeoPKPD.DosingRegimen(NeoPKPD.QD(), dose, 28)

    return NeoPKPD.TreatmentArm(name, model_spec, regimen; n_subjects=n_subjects)
end

function _parse_endpoint(endpoint_dict)
    name = Symbol(endpoint_dict["name"])
    endpoint_type = lowercase(get(endpoint_dict, "type", "pk"))

    if endpoint_type == "pk"
        metric = Symbol(get(endpoint_dict, "metric", "auc_0_inf"))
        return NeoPKPD.PKEndpoint(name; metric=metric)
    elseif endpoint_type == "pd"
        metric = Symbol(get(endpoint_dict, "metric", "change_from_baseline"))
        return NeoPKPD.PDEndpoint(name; metric=metric)
    elseif endpoint_type == "safety"
        return NeoPKPD.SafetyEndpoint(name)
    else
        return NeoPKPD.PKEndpoint(name)
    end
end

# ============================================================================
# Read CDISC Command
# ============================================================================

const READ_CDISC_HELP = """
Read CDISC/SDTM data and convert to NeoPKPD format.

Arguments:
  --pc PATH           Path to PC (Pharmacokinetic Concentrations) domain CSV (required)
  --ex PATH           Path to EX (Exposure) domain CSV (optional)
  --dm PATH           Path to DM (Demographics) domain CSV (optional)
  --out PATH          Output path for converted data (required)

Examples:
  neopkpd read-cdisc --pc pc.csv --out data.json
  neopkpd read-cdisc --pc pc.csv --ex ex.csv --dm dm.csv --out data.json
"""

function cmd_read_cdisc(args)
    pc_path = args["pc"]
    ex_path = get(args, "ex", nothing)
    dm_path = get(args, "dm", nothing)
    out = args["out"]

    if !isfile(pc_path)
        _die("PC file not found: $pc_path")
    end

    _info("Reading CDISC data...")

    # Read PC domain
    pc_data = NeoPKPD.read_cdisc_pc(pc_path)

    # Read optional domains
    ex_data = ex_path !== nothing && isfile(ex_path) ? NeoPKPD.read_cdisc_ex(ex_path) : nothing
    dm_data = dm_path !== nothing && isfile(dm_path) ? NeoPKPD.read_cdisc_dm(dm_path) : nothing

    # Convert to observed data format
    observed = NeoPKPD.cdisc_to_observed_data(pc_data, ex_data, dm_data)

    # Build output
    output = Dict(
        "n_subjects" => NeoPKPD.n_subjects(observed),
        "n_observations" => NeoPKPD.n_observations(observed),
        "subjects" => [
            Dict(
                "subject_id" => s.subject_id,
                "n_observations" => length(s.times),
                "time_range" => [minimum(s.times), maximum(s.times)]
            ) for s in observed.subjects
        ]
    )

    if dm_data !== nothing
        output["demographics_available"] = true
    end

    _write_json(out, output)
    _info("CDISC data converted and written to: $out")
    _info("Subjects: $(NeoPKPD.n_subjects(observed)), Observations: $(NeoPKPD.n_observations(observed))")
end

# ============================================================================
# Help Command
# ============================================================================

const MAIN_HELP = """
NeoPKPD - Professional PK/PD Simulation Platform

USAGE:
    neopkpd <command> [options]

COMMANDS:
    version          Display version information
    simulate         Run a PK or PKPD simulation from a JSON spec
    population       Run a population simulation with IIV/IOV
    sensitivity      Run parameter sensitivity analysis
    metrics          Compute PK/PD metrics from simulation results
    replay           Replay a simulation from an artifact file
    validate-golden  Run golden validation tests
    import           Import NONMEM/Monolix model files
    nca              Run Non-Compartmental Analysis
    vpc              Run Visual Predictive Check
    estimate         Run parameter estimation (NLME)
    read-cdisc       Read CDISC/SDTM data files
    trial            Run clinical trial simulation
    help             Show this help message

SUPPORTED PK MODELS:
    OneCompIVBolus, OneCompOralFirstOrder, TwoCompIVBolus, TwoCompOral,
    ThreeCompIVBolus, TransitAbsorption, MichaelisMentenElimination

SUPPORTED PD MODELS:
    DirectEmax, SigmoidEmax, IndirectResponseTurnover, BiophaseEquilibration

Use 'neopkpd help <command>' for detailed help on each command.

EXAMPLES:
    neopkpd version
    neopkpd simulate --spec pk_spec.json --out result.json
    neopkpd population --spec pop_spec.json --out pop_result.json
    neopkpd metrics --artifact result.json --metric cmax
    neopkpd replay --artifact validation/golden/pk_iv_bolus.json
    neopkpd import --file model.ctl --out spec.json
    neopkpd nca --data pk_data.csv --dose 100 --out nca.json
    neopkpd vpc --observed data.csv --pop-spec pop.json --out vpc.json
    neopkpd estimate --data data.csv --model-spec model.json --out est.json

For more information, visit: https://github.com/neopkpd/neopkpd
"""

function cmd_help(command::Union{String,Nothing})
    if command === nothing
        println(MAIN_HELP)
    elseif command == "version"
        println(VERSION_HELP)
    elseif command == "simulate"
        println(SIMULATE_HELP)
    elseif command == "population"
        println(POPULATION_HELP)
    elseif command == "sensitivity"
        println(SENSITIVITY_HELP)
    elseif command == "metrics"
        println(METRICS_HELP)
    elseif command == "replay"
        println(REPLAY_HELP)
    elseif command == "validate-golden"
        println(VALIDATE_GOLDEN_HELP)
    elseif command == "import"
        println(IMPORT_HELP)
    elseif command == "nca"
        println(NCA_HELP)
    elseif command == "vpc"
        println(VPC_HELP)
    elseif command == "estimate"
        println(ESTIMATE_HELP)
    elseif command == "read-cdisc"
        println(READ_CDISC_HELP)
    elseif command == "trial"
        println(TRIAL_HELP)
    else
        println("Unknown command: $command")
        println()
        println(MAIN_HELP)
    end
end

# ============================================================================
# Main Entry Point
# ============================================================================

function main()
    if length(ARGS) < 1
        println(MAIN_HELP)
        return
    end

    cmd = ARGS[1]
    rest = ARGS[2:end]

    if cmd == "version"
        cmd_version()

    elseif cmd == "help"
        if length(rest) > 0
            cmd_help(rest[1])
        else
            cmd_help(nothing)
        end

    elseif cmd == "replay"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--artifact"
                required = true
                help = "Path to the artifact JSON file"
            "--out"
                required = false
                help = "Output path for replayed artifact"
        end
        args = parse_args(rest, rs)
        cmd_replay(args)

    elseif cmd == "validate-golden"
        cmd_validate_golden()

    elseif cmd == "simulate"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to simulation specification JSON"
            "--out"
                required = true
                help = "Output path for results"
            "--format"
                required = false
                default = "artifact"
                help = "Output format: artifact or simple"
        end
        args = parse_args(rest, rs)
        cmd_simulate(args)

    elseif cmd == "population"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to population specification JSON"
            "--out"
                required = true
                help = "Output path for results"
            "--format"
                required = false
                default = "artifact"
                help = "Output format: artifact or simple"
        end
        args = parse_args(rest, rs)
        cmd_population(args)

    elseif cmd == "sensitivity"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to sensitivity specification JSON"
            "--out"
                required = true
                help = "Output path for results"
        end
        args = parse_args(rest, rs)
        cmd_sensitivity(args)

    elseif cmd == "metrics"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--artifact"
                required = true
                help = "Path to artifact or result JSON"
            "--observation"
                required = false
                default = "conc"
                help = "Observation to analyze"
            "--metric"
                required = true
                help = "Metric: cmax, auc, emin, time_below, auc_above_baseline"
            "--threshold"
                required = false
                help = "Threshold for time_below or auc_above_baseline"
        end
        args = parse_args(rest, rs)
        cmd_metrics(args)

    elseif cmd == "import"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--file"
                required = true
                help = "Path to .ctl (NONMEM) or .mlxtran (Monolix) file"
            "--format"
                required = false
                help = "Input format: nonmem or monolix (auto-detected if not specified)"
            "--out"
                required = true
                help = "Output path for NeoPKPD JSON spec"
        end
        args = parse_args(rest, rs)
        cmd_import(args)

    elseif cmd == "nca"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--data"
                required = true
                help = "Path to data file (CSV with time, conc columns)"
            "--dose"
                required = true
                help = "Dose amount for CL/F and Vz/F calculation"
            "--out"
                required = true
                help = "Output path for NCA results"
            "--method"
                required = false
                default = "linear_up_log_down"
                help = "AUC method: linear, log_linear, linear_up_log_down"
            "--tau"
                required = false
                help = "Dosing interval for steady-state metrics"
        end
        args = parse_args(rest, rs)
        cmd_nca(args)

    elseif cmd == "vpc"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--observed"
                required = true
                help = "Path to observed data CSV"
            "--pop-spec"
                required = true
                help = "Path to population specification JSON"
            "--out"
                required = true
                help = "Output path for VPC results"
            "--n-sim"
                required = false
                default = "200"
                help = "Number of simulations"
            "--n-bins"
                required = false
                default = "8"
                help = "Number of time bins"
            "--pi"
                required = false
                default = "0.05,0.50,0.95"
                help = "Prediction interval levels, comma-separated"
            "--pcvpc"
                action = :store_true
                help = "Use prediction-corrected VPC"
        end
        args = parse_args(rest, rs)
        cmd_vpc(args)

    elseif cmd == "estimate"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--data"
                required = true
                help = "Path to observed data CSV"
            "--model-spec"
                required = true
                help = "Path to model specification JSON"
            "--out"
                required = true
                help = "Output path for estimation results"
            "--method"
                required = false
                default = "foce"
                help = "Estimation method: laplacian, foce, saem"
            "--max-iter"
                required = false
                default = "500"
                help = "Maximum iterations"
            "--compute-se"
                required = false
                default = "true"
                help = "Compute standard errors"
        end
        args = parse_args(rest, rs)
        cmd_estimate(args)

    elseif cmd == "read-cdisc"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--pc"
                required = true
                help = "Path to PC (Pharmacokinetic Concentrations) domain CSV"
            "--ex"
                required = false
                help = "Path to EX (Exposure) domain CSV"
            "--dm"
                required = false
                help = "Path to DM (Demographics) domain CSV"
            "--out"
                required = true
                help = "Output path for converted data"
        end
        args = parse_args(rest, rs)
        cmd_read_cdisc(args)

    elseif cmd == "trial"
        rs = ArgParseSettings()
        @add_arg_table rs begin
            "--spec"
                required = true
                help = "Path to trial specification JSON"
            "--out"
                required = true
                help = "Output path for trial results"
            "--replicates"
                required = false
                default = "1"
                help = "Number of trial replicates for power analysis"
        end
        args = parse_args(rest, rs)
        cmd_trial(args)

    else
        _die("Unknown command: $cmd. Use 'neopkpd help' for usage.")
    end
end

# Run when invoked as script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
