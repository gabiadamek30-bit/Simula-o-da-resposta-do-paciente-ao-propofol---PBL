using Pkg
Pkg.activate("packages/core")
Pkg.instantiate()

using NeoPKPD
using JSON

function write(path, artifact)
    open(path, "w") do io
        JSON.print(io, artifact)
    end
end

function gen_pk_iv_bolus()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [
            DoseEvent(0.0, 60.0),
            DoseEvent(0.0, 40.0),   # duplicate time to lock semantics summing
            DoseEvent(12.0, 25.0),
        ],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_pk_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_then_pd_direct_emax()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_then_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        DirectEmax(),
        "golden_emax",
        DirectEmaxParams(10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_coupled_pkpd_turnover_oral()
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "golden_coupled_oral",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_population_pk_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(7777), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

function gen_sensitivity_single_iv()
    spec = ModelSpec(
        OneCompIVBolus(),
        "golden_sens_single",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan("CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)])

    res = run_sensitivity(spec, grid, solver; plan = plan, observation = :conc)

    return serialize_sensitivity_execution(model_spec = spec, grid = grid, solver = solver, result = res)
end

function gen_sensitivity_population_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_sens_pop_base",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    iiv = IIVSpec(LogNormalIIV(), Dict(:CL => 0.2, :V => 0.1), UInt64(7777), 5)
    pop = PopulationSpec(base, iiv, nothing, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    plan = PerturbationPlan("CL_up_10pct", [Perturbation(RelativePerturbation(), :CL, 0.10)])

    res = run_population_sensitivity(pop, grid, solver; plan = plan, observation = :conc, probs = [0.05, 0.95])

    return serialize_population_sensitivity_execution(population_spec = pop, grid = grid, solver = solver, result = res)
end

function gen_population_iov_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_iov_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    iov = IOVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(1234), OccasionDefinition(:dose_times))
    pop = PopulationSpec(base, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

function gen_population_iov_pkpd()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_iov_pkpd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    iov = IOVSpec(LogNormalIIV(), Dict(:CL => 0.3), UInt64(2222), OccasionDefinition(:dose_times))
    pop = PopulationSpec(pk, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver; pd_spec=pd)

    return serialize_population_execution(
        population_spec=pop,
        grid=grid,
        solver=solver,
        result=res,
        pd_spec=pd,
    )
end


function gen_population_time_varying_covariate_iv()
    base = ModelSpec(
        OneCompIVBolus(),
        "golden_tv_cov_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    cm = CovariateModel(
        "cl_tv",
        [CovariateEffect(LinearCovariate(), :CL, :CLMULT, 1.0, 1.0)],
    )

    tv = TimeVaryingCovariates(Dict(
        :CLMULT => TimeCovariateSeries(StepTimeCovariate(), [0.0, 10.0], [1.0, 2.0]),
    ))

    covs = [IndividualCovariates(Dict{Symbol,Float64}(), tv)]
    pop = PopulationSpec(base, nothing, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    return serialize_population_execution(population_spec=pop, grid=grid, solver=solver, result=res)
end

# =====================================================
# New PK Models
# =====================================================

function gen_pk_twocomp_iv()
    pk = ModelSpec(
        TwoCompIVBolus(),
        "golden_pk_twocomp_iv",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),  # CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0), DoseEvent(12.0, 300.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_twocomp_oral()
    pk = ModelSpec(
        TwoCompOral(),
        "golden_pk_twocomp_oral",
        TwoCompOralParams(1.0, 10.0, 50.0, 5.0, 100.0),  # Ka=1, CL=10, V1=50, Q=5, V2=100
        [DoseEvent(0.0, 500.0), DoseEvent(12.0, 300.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_threecomp_iv()
    pk = ModelSpec(
        ThreeCompIVBolus(),
        "golden_pk_threecomp_iv",
        ThreeCompIVBolusParams(10.0, 50.0, 10.0, 80.0, 2.0, 200.0),
        [DoseEvent(0.0, 1000.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_transit_absorption()
    pk = ModelSpec(
        TransitAbsorption(),
        "golden_pk_transit",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),  # N=5, Ktr=2, Ka=1, CL=10, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

function gen_pk_michaelis_menten()
    pk = ModelSpec(
        MichaelisMentenElimination(),
        "golden_pk_mm",
        MichaelisMentenEliminationParams(100.0, 5.0, 50.0),  # Vmax=100, Km=5, V=50
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(pk, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res)
end

# =====================================================
# New PD Models
# =====================================================

function gen_pkpd_sigmoid_emax()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_for_semax",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        SigmoidEmax(),
        "golden_semax",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 2.0),  # E0=10, Emax=40, EC50=0.8, gamma=2
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_pkpd_biophase_equilibration()
    pk = ModelSpec(
        OneCompIVBolus(),
        "golden_pk_for_biophase",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    pd = PDSpec(
        BiophaseEquilibration(),
        "golden_biophase",
        BiophaseEquilibrationParams(0.5, 10.0, 40.0, 0.8),  # ke0=0.5, E0=10, Emax=40, EC50=0.8
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

# =====================================================
# Integration: New PK + PD combinations
# =====================================================

function gen_twocomp_with_sigmoid_emax()
    pk = ModelSpec(
        TwoCompIVBolus(),
        "golden_twocomp_semax",
        TwoCompIVBolusParams(10.0, 50.0, 5.0, 100.0),
        [DoseEvent(0.0, 500.0)],
    )

    pd = PDSpec(
        SigmoidEmax(),
        "golden_semax_2c",
        SigmoidEmaxParams(0.0, 100.0, 2.0, 2.0),
        :conc,
        :effect,
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

function gen_transit_with_turnover()
    pk = ModelSpec(
        TransitAbsorption(),
        "golden_transit_turnover",
        TransitAbsorptionParams(5, 2.0, 1.0, 10.0, 50.0),
        [DoseEvent(0.0, 500.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "golden_turnover_transit",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 1.0),
        :conc,
        :response,
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    return serialize_execution(model_spec=pk, grid=grid, solver=solver, result=res, pd_spec=pd)
end

# =====================================================
# TMDD (Target-Mediated Drug Disposition) Models
# Industry-standard models for monoclonal antibodies
# =====================================================

function gen_tmdd_qss_iv()
    # Quasi-Steady-State 2-compartment TMDD - industry standard for mAbs
    # Parameters typical for IgG monoclonal antibody
    spec = TMDDSpec(
        TwoCptTMDD(QSS, IVBolus),
        "golden_tmdd_qss_iv",
        TwoCptTMDDParams(
            0.2,     # CL = 0.2 L/day (typical for IgG)
            3.0,     # V1 = 3 L (plasma volume)
            2.5,     # V2 = 2.5 L (peripheral)
            0.5,     # Q = 0.5 L/day
            1.0,     # KSS = 1.0 nM
            0.1,     # kint = 0.1 1/day
            0.1,     # ksyn = 0.1 nM/day
            0.1,     # kdeg = 0.1 1/day
            1.0,     # R0 = 1.0 nM
        ),
        [DoseEvent(0.0, 200.0)],  # 200 mg IV bolus
    )

    grid = SimGrid(0.0, 28.0, collect(0.0:0.5:28.0))  # 28 days
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_full_iv()
    # Full mechanistic TMDD (Mager-Jusko model)
    spec = TMDDSpec(
        TwoCptTMDD(FullTMDD, IVBolus),
        "golden_tmdd_full_iv",
        TwoCptTMDDParams(
            0.25,    # CL = 0.25 L/day
            3.5,     # V1 = 3.5 L
            3.0,     # V2 = 3.0 L
            0.6,     # Q = 0.6 L/day
            0.5,     # KSS = 0.5 nM
            0.15,    # kint = 0.15 1/day
            0.15,    # ksyn = 0.15 nM/day
            0.15,    # kdeg = 0.15 1/day
            1.0,     # R0 = 1.0 nM
        ),
        [DoseEvent(0.0, 150.0)],  # 150 mg IV bolus
    )

    grid = SimGrid(0.0, 21.0, collect(0.0:0.5:21.0))  # 21 days
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_qss_sc()
    # QSS TMDD with subcutaneous administration
    spec = TMDDSpec(
        TwoCptTMDD(QSS, Subcutaneous),
        "golden_tmdd_qss_sc",
        TwoCptTMDDParams(
            0.2,     # CL
            3.0,     # V1
            2.5,     # V2
            0.5,     # Q
            1.0,     # KSS
            0.1,     # kint
            0.1,     # ksyn
            0.1,     # kdeg
            1.0,     # R0
            0.3,     # ka = 0.3 1/day (SC absorption)
            0.7,     # F = 0.7 (70% bioavailability)
            0.5,     # Tlag = 0.5 days
        ),
        [DoseEvent(0.0, 300.0)],  # 300 mg SC
    )

    grid = SimGrid(0.0, 28.0, collect(0.0:0.5:28.0))
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_rapid_binding()
    # Rapid binding (Wagner) approximation
    spec = TMDDSpec(
        TwoCptTMDD(RapidBinding, IVBolus),
        "golden_tmdd_rapid_binding",
        TwoCptTMDDParams(
            0.3,     # CL
            4.0,     # V1
            3.0,     # V2
            0.8,     # Q
            2.0,     # KSS
            0.2,     # kint
            0.2,     # ksyn
            0.2,     # kdeg
            1.0,     # R0
        ),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 14.0, collect(0.0:0.25:14.0))  # 14 days
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_fcrn()
    # FcRn recycling model for IgG half-life
    spec = TMDDSpec(
        TwoCptTMDDFcRn(QSS, IVBolus),
        "golden_tmdd_fcrn",
        TwoCptTMDDFcRnParams(
            3.0,     # V1
            2.5,     # V2
            0.5,     # Q
            0.3,     # CLup (pinocytic uptake)
            0.9,     # FR (90% FcRn recycling)
            1.0,     # KSS
            0.1,     # kint
            0.1,     # ksyn
            0.1,     # kdeg
            1.0,     # R0
        ),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 42.0, collect(0.0:1.0:42.0))  # 42 days (6 weeks)
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_multiple_dose()
    # Multiple dose TMDD - Q3W dosing typical for mAbs
    spec = TMDDSpec(
        TwoCptTMDD(QSS, IVBolus),
        "golden_tmdd_multiple_dose",
        TwoCptTMDDParams(
            0.2,     # CL
            3.0,     # V1
            2.5,     # V2
            0.5,     # Q
            1.0,     # KSS
            0.1,     # kint
            0.1,     # ksyn
            0.1,     # kdeg
            1.0,     # R0
        ),
        [
            DoseEvent(0.0, 200.0),    # Dose 1: Day 0
            DoseEvent(21.0, 200.0),   # Dose 2: Day 21 (Q3W)
            DoseEvent(42.0, 200.0),   # Dose 3: Day 42
            DoseEvent(63.0, 200.0),   # Dose 4: Day 63
        ],
    )

    grid = SimGrid(0.0, 84.0, collect(0.0:1.0:84.0))  # 12 weeks
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_onecomp_qss()
    # One-compartment QSS TMDD for simpler biologics
    spec = TMDDSpec(
        OneCptTMDD(QSS, IVBolus),
        "golden_tmdd_onecomp_qss",
        OneCptTMDDParams(
            0.3,     # CL = 0.3 L/day
            5.0,     # V = 5 L
            0.5,     # KSS = 0.5 nM
            0.2,     # kint = 0.2 1/day
            0.2,     # ksyn = 0.2 nM/day
            0.2,     # kdeg = 0.2 1/day
            1.0,     # R0 = 1.0 nM
        ),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 14.0, collect(0.0:0.25:14.0))
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function gen_tmdd_high_target()
    # TMDD with high baseline target (saturable kinetics more prominent)
    spec = TMDDSpec(
        TwoCptTMDD(QSS, IVBolus),
        "golden_tmdd_high_target",
        TwoCptTMDDParams(
            0.15,    # CL (lower linear clearance)
            3.0,     # V1
            2.5,     # V2
            0.4,     # Q
            0.1,     # KSS = 0.1 nM (high affinity)
            0.3,     # kint (faster internalization)
            1.0,     # ksyn = 1.0 nM/day (high synthesis)
            0.1,     # kdeg = 0.1 1/day
            10.0,    # R0 = 10.0 nM (high baseline target)
        ),
        [DoseEvent(0.0, 50.0)],  # Lower dose to show nonlinear kinetics
    )

    grid = SimGrid(0.0, 21.0, collect(0.0:0.5:21.0))
    solver = SolverSpec(:Rodas5, 1e-10, 1e-12, 10^7)

    res = solve_tmdd(spec, grid, solver)

    return serialize_tmdd_execution(tmdd_spec=spec, grid=grid, solver=solver, result=res)
end

function main()
    mkpath("validation/golden")

    artifacts = Dict(
        # Original goldens
        "pk_iv_bolus.json" => gen_pk_iv_bolus(),
        "pk_oral.json" => gen_pk_oral(),
        "pk_then_pd_direct_emax.json" => gen_pk_then_pd_direct_emax(),
        "pkpd_coupled_turnover_oral.json" => gen_coupled_pkpd_turnover_oral(),
        "population_pk_iv.json" => gen_population_pk_iv(),
        "sensitivity_single_iv.json" => gen_sensitivity_single_iv(),
        "sensitivity_population_iv.json" => gen_sensitivity_population_iv(),
        "population_iov_iv.json" => gen_population_iov_iv(),
        "population_iov_pkpd.json" => gen_population_iov_pkpd(),
        "population_time_varying_cov_iv.json" => gen_population_time_varying_covariate_iv(),

        # New PK model goldens
        "pk_twocomp_iv.json" => gen_pk_twocomp_iv(),
        "pk_twocomp_oral.json" => gen_pk_twocomp_oral(),
        "pk_threecomp_iv.json" => gen_pk_threecomp_iv(),
        "pk_transit_absorption.json" => gen_pk_transit_absorption(),
        "pk_michaelis_menten.json" => gen_pk_michaelis_menten(),

        # New PD model goldens
        "pkpd_sigmoid_emax.json" => gen_pkpd_sigmoid_emax(),
        "pkpd_biophase_equilibration.json" => gen_pkpd_biophase_equilibration(),

        # Integration goldens
        "pkpd_twocomp_sigmoid_emax.json" => gen_twocomp_with_sigmoid_emax(),
        "pkpd_transit_turnover.json" => gen_transit_with_turnover(),

        # TMDD (Target-Mediated Drug Disposition) goldens
        "tmdd_qss_iv.json" => gen_tmdd_qss_iv(),
        "tmdd_full_iv.json" => gen_tmdd_full_iv(),
        "tmdd_qss_sc.json" => gen_tmdd_qss_sc(),
        "tmdd_rapid_binding.json" => gen_tmdd_rapid_binding(),
        "tmdd_fcrn.json" => gen_tmdd_fcrn(),
        "tmdd_multiple_dose.json" => gen_tmdd_multiple_dose(),
        "tmdd_onecomp_qss.json" => gen_tmdd_onecomp_qss(),
        "tmdd_high_target.json" => gen_tmdd_high_target(),
    )

    for (fname, art) in artifacts
        path = joinpath("validation/golden", fname)
        write(path, art)
        println("Wrote: $(path)")
    end
end

main()
