# Coupled PK-PD model tests

@testset "PKPD coupling with DirectEmax using IV bolus analytic reference" begin
    pk = ModelSpec(
        OneCompIVBolus(), "pk_iv", OneCompIVBolusParams(5.0, 50.0), [DoseEvent(0.0, 100.0)]
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(DirectEmax(), "pd_emax", DirectEmaxParams(10.0, 40.0, 0.8), :conc, :effect)

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    CL = pk.params.CL
    V = pk.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        e_ref = direct_emax(c_ref, pd.params.E0, pd.params.Emax, pd.params.EC50)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10, atol=1e-12)
    end
end

@testset "IndirectResponseTurnover coupled: Imax=0 matches analytic turnover and PK analytic" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_for_pd",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax, no drug effect
            1.0,   # IC50, irrelevant here
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in pk.doses)
            continue  # skip discontinuities (left-continuous solver output)
        end

        c_ref = analytic_onecomp_ivbolus_conc(t, pk.doses, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "IndirectResponseTurnover coupled: inhibition raises response above baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    # Baseline at steady state to make interpretation clean
    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_with_effect",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]
    @test maximum(r) > Rss
end

@testset "Coupled engine generalization: oral PK with Imax=0 matches analytic PK and analytic turnover" begin
    pk = ModelSpec(
        OneCompOralFirstOrder(),
        "pk_oral_for_pd",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        IndirectResponseTurnover(),
        "turnover_no_effect_oral",
        IndirectResponseTurnoverParams(
            10.0,  # Kin
            0.5,   # Kout
            15.0,  # R0
            0.0,   # Imax
            1.0,   # IC50
        ),
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    Ka = pk.params.Ka
    CL = pk.params.CL
    V = pk.params.V

    Kin = pd.params.Kin
    Kout = pd.params.Kout
    R0 = pd.params.R0

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, pk.doses, Ka, CL, V)
        r_ref = analytic_turnover_R(t, Kin, Kout, R0)

        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
        @test isapprox(res.observations[:response][i], r_ref; rtol=1e-8, atol=1e-10)
    end
end

# =============================================================================
# IRM-I Tests: Inhibition of Kin (Production)
# dR/dt = Kin × (1 - I(C)) - Kout × R
# =============================================================================

@testset "IRM-I coupled: Imax=0 matches analytic turnover (no drug effect)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm1",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    R0 = Kin / Kout  # At steady state

    pd = PDSpec(
        IndirectResponseIRM1(),
        "irm1_no_effect",
        IndirectResponseIRM1Params(Kin, Kout, R0, 0.0, 1.0),  # Imax=0
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    # With Imax=0, response should stay at baseline (R0)
    for (i, t) in enumerate(res.t)
        if t == 0.0
            continue  # Skip initial time
        end
        @test isapprox(res.observations[:response][i], R0; rtol=1e-6)
    end
end

@testset "IRM-I coupled: Inhibition of Kin decreases response below baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm1_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseIRM1(),
        "irm1_with_effect",
        IndirectResponseIRM1Params(Kin, Kout, Rss, 0.8, 0.5),  # Imax=0.8
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]

    # IRM-I inhibits production → response should decrease below baseline
    @test minimum(r) < Rss

    # Response should return toward baseline as drug washes out
    @test r[end] > minimum(r)
end

# =============================================================================
# IRM-II Tests: Stimulation of Kin (Production)
# dR/dt = Kin × (1 + S(C)) - Kout × R
# =============================================================================

@testset "IRM-II coupled: Smax=0 matches analytic turnover (no drug effect)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm2",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    R0 = Kin / Kout

    pd = PDSpec(
        IndirectResponseIRM2(),
        "irm2_no_effect",
        IndirectResponseIRM2Params(Kin, Kout, R0, 0.0, 1.0),  # Smax=0
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    # With Smax=0, response should stay at baseline (R0)
    for (i, t) in enumerate(res.t)
        if t == 0.0
            continue
        end
        @test isapprox(res.observations[:response][i], R0; rtol=1e-6)
    end
end

@testset "IRM-II coupled: Stimulation of Kin increases response above baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm2_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseIRM2(),
        "irm2_with_effect",
        IndirectResponseIRM2Params(Kin, Kout, Rss, 2.0, 0.5),  # Smax=2.0 (200% stimulation)
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]

    # IRM-II stimulates production → response should increase above baseline
    @test maximum(r) > Rss

    # Response should return toward baseline as drug washes out
    @test r[end] < maximum(r)
end

# =============================================================================
# IRM-IV Tests: Stimulation of Kout (Elimination)
# dR/dt = Kin - Kout × (1 + S(C)) × R
# =============================================================================

@testset "IRM-IV coupled: Smax=0 matches analytic turnover (no drug effect)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm4",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    R0 = Kin / Kout

    pd = PDSpec(
        IndirectResponseIRM4(),
        "irm4_no_effect",
        IndirectResponseIRM4Params(Kin, Kout, R0, 0.0, 1.0),  # Smax=0
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :response)

    # With Smax=0, response should stay at baseline (R0)
    for (i, t) in enumerate(res.t)
        if t == 0.0
            continue
        end
        @test isapprox(res.observations[:response][i], R0; rtol=1e-6)
    end
end

@testset "IRM-IV coupled: Stimulation of Kout decreases response below baseline" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv_irm4_effect",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseIRM4(),
        "irm4_with_effect",
        IndirectResponseIRM4Params(Kin, Kout, Rss, 2.0, 0.5),  # Smax=2.0
        :conc,
        :response,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    r = res.observations[:response]

    # IRM-IV stimulates elimination → response should decrease below baseline
    @test minimum(r) < Rss

    # Response should return toward baseline as drug washes out
    @test r[end] > minimum(r)
end

# =============================================================================
# IRM Type Comparison Tests
# =============================================================================

@testset "IRM-III (IndirectResponseTurnover) and IRM-III alias produce identical results" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    Kin = 10.0
    Kout = 0.5
    R0 = 15.0

    # Using original type name
    pd_original = PDSpec(
        IndirectResponseTurnover(),
        "irm3_original",
        IndirectResponseTurnoverParams(Kin, Kout, R0, 0.5, 1.0),
        :conc,
        :response,
    )

    # Using IRM3 alias (should be identical)
    pd_alias = PDSpec(
        IndirectResponseIRM3(),
        "irm3_alias",
        IndirectResponseIRM3Params(Kin, Kout, R0, 0.5, 1.0),
        :conc,
        :response,
    )

    res_original = simulate_pkpd_coupled(pk, pd_original, grid, solver)
    res_alias = simulate_pkpd_coupled(pk, pd_alias, grid, solver)

    @test res_original.t == res_alias.t
    @test res_original.observations[:response] ≈ res_alias.observations[:response]
end

@testset "IRM-I and IRM-III have opposite effects on response direction" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout
    Imax = 0.8
    IC50 = 0.5

    # IRM-I: Inhibition of Kin → Response decreases
    pd_irm1 = PDSpec(
        IndirectResponseIRM1(),
        "irm1",
        IndirectResponseIRM1Params(Kin, Kout, Rss, Imax, IC50),
        :conc,
        :response,
    )

    # IRM-III: Inhibition of Kout → Response increases
    pd_irm3 = PDSpec(
        IndirectResponseTurnover(),
        "irm3",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, Imax, IC50),
        :conc,
        :response,
    )

    res_irm1 = simulate_pkpd_coupled(pk, pd_irm1, grid, solver)
    res_irm3 = simulate_pkpd_coupled(pk, pd_irm3, grid, solver)

    r1 = res_irm1.observations[:response]
    r3 = res_irm3.observations[:response]

    # IRM-I should decrease below baseline
    @test minimum(r1) < Rss

    # IRM-III should increase above baseline
    @test maximum(r3) > Rss
end

@testset "IRM-II and IRM-IV have opposite effects on response direction" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout
    Smax = 2.0
    SC50 = 0.5

    # IRM-II: Stimulation of Kin → Response increases
    pd_irm2 = PDSpec(
        IndirectResponseIRM2(),
        "irm2",
        IndirectResponseIRM2Params(Kin, Kout, Rss, Smax, SC50),
        :conc,
        :response,
    )

    # IRM-IV: Stimulation of Kout → Response decreases
    pd_irm4 = PDSpec(
        IndirectResponseIRM4(),
        "irm4",
        IndirectResponseIRM4Params(Kin, Kout, Rss, Smax, SC50),
        :conc,
        :response,
    )

    res_irm2 = simulate_pkpd_coupled(pk, pd_irm2, grid, solver)
    res_irm4 = simulate_pkpd_coupled(pk, pd_irm4, grid, solver)

    r2 = res_irm2.observations[:response]
    r4 = res_irm4.observations[:response]

    # IRM-II should increase above baseline
    @test maximum(r2) > Rss

    # IRM-IV should decrease below baseline
    @test minimum(r4) < Rss
end

# =============================================================================
# Transit Compartment PD Tests
# Signal transduction delay model
# =============================================================================

@testset "Transit compartment PD: baseline with no drug effect (Emax=0)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    E0 = 10.0  # Baseline
    N = 3
    ktr = 0.5

    pd = PDSpec(
        TransitCompartmentPD(),
        "transit_no_effect",
        TransitCompartmentPDParams(N, ktr, E0, 0.0, 1.0, 1.0),  # Emax=0
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test haskey(res.states, :A1)
    @test haskey(res.states, :A2)
    @test haskey(res.states, :A3)

    # With Emax=0, effect should remain at baseline E0
    for (i, t) in enumerate(res.t)
        @test isapprox(res.observations[:effect][i], E0; rtol=1e-4)
    end
end

@testset "Transit compartment PD: delayed response characteristics" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:0.5:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    E0 = 10.0
    Emax = 20.0
    EC50 = 1.0
    gamma = 1.0
    N = 5
    ktr = 0.3

    pd = PDSpec(
        TransitCompartmentPD(),
        "transit_delayed",
        TransitCompartmentPDParams(N, ktr, E0, Emax, EC50, gamma),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    effect = res.observations[:effect]
    conc = res.observations[:conc]

    # Find time of peak concentration
    t_max_conc = res.t[argmax(conc)]

    # Find time of peak effect
    t_max_effect = res.t[argmax(effect)]

    # Effect should peak AFTER concentration (delayed response)
    @test t_max_effect > t_max_conc

    # Effect should increase above baseline
    @test maximum(effect) > E0

    # Check MTT is approximately correct
    # MTT = (N + 1) / ktr
    expected_mtt = (N + 1) / ktr
    @test isapprox(res.metadata["pd_params"]["MTT"], expected_mtt; rtol=1e-10)
end

@testset "Transit compartment PD: increasing N increases delay" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:0.5:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    E0 = 10.0
    Emax = 15.0
    EC50 = 1.0
    gamma = 1.0
    ktr = 0.3

    # N=2 (fewer compartments, less delay)
    pd_n2 = PDSpec(
        TransitCompartmentPD(),
        "transit_n2",
        TransitCompartmentPDParams(2, ktr, E0, Emax, EC50, gamma),
        :conc,
        :effect,
    )

    # N=6 (more compartments, more delay)
    pd_n6 = PDSpec(
        TransitCompartmentPD(),
        "transit_n6",
        TransitCompartmentPDParams(6, ktr, E0, Emax, EC50, gamma),
        :conc,
        :effect,
    )

    res_n2 = simulate_pkpd_coupled(pk, pd_n2, grid, solver)
    res_n6 = simulate_pkpd_coupled(pk, pd_n6, grid, solver)

    effect_n2 = res_n2.observations[:effect]
    effect_n6 = res_n6.observations[:effect]

    # Find time of peak effect for each
    t_peak_n2 = res_n2.t[argmax(effect_n2)]
    t_peak_n6 = res_n6.t[argmax(effect_n6)]

    # More compartments → longer delay to peak
    @test t_peak_n6 > t_peak_n2

    # MTT should scale correctly
    mtt_n2 = (2 + 1) / ktr
    mtt_n6 = (6 + 1) / ktr
    @test mtt_n6 > mtt_n2
end

@testset "Transit compartment PD: single compartment (N=1)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        TransitCompartmentPD(),
        "transit_n1",
        TransitCompartmentPDParams(1, 0.5, 10.0, 10.0, 1.0, 1.0),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.states, :A1)
    @test !haskey(res.states, :A2)  # Only one transit compartment

    # Effect should equal A1 for N=1
    @test res.observations[:effect] ≈ res.states[:A1]
end

# =============================================================================
# Disease Progression PD Tests
# Tumor growth dynamics with drug effect
# =============================================================================

@testset "Disease progression: Exponential growth without drug (kdrug=0)" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    S0 = 10.0
    kgrow = 0.05  # 5% growth per time unit

    pd = PDSpec(
        DiseaseProgressionPD(ExponentialGrowth),
        "tumor_exp",
        DiseaseProgressionPDParams(S0, kgrow, 100.0, 0.0, 0.0),  # kdrug=0
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :tumor)
    @test haskey(res.states, :S)

    # Without drug, tumor should grow exponentially: S(t) = S0 * exp(kgrow * t)
    for (i, t) in enumerate(res.t)
        expected_S = S0 * exp(kgrow * t)
        @test isapprox(res.observations[:tumor][i], expected_S; rtol=1e-4)
    end
end

@testset "Disease progression: Gompertz growth without drug" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 100.0, collect(0.0:2.0:100.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    S0 = 10.0
    kgrow = 0.05
    Smax = 100.0

    pd = PDSpec(
        DiseaseProgressionPD(GompertzGrowth),
        "tumor_gompertz",
        DiseaseProgressionPDParams(S0, kgrow, Smax, 0.0, 0.0),  # kdrug=0
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    tumor = res.observations[:tumor]

    # Gompertz should approach Smax asymptotically
    @test tumor[end] > S0
    @test tumor[end] < Smax

    # Growth should slow as tumor approaches Smax
    growth_early = tumor[5] - tumor[1]
    growth_late = tumor[end] - tumor[end-4]
    @test growth_late < growth_early  # Growth slows over time
end

@testset "Disease progression: Logistic growth without drug" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 100.0, collect(0.0:2.0:100.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    S0 = 10.0
    kgrow = 0.1
    Smax = 100.0

    pd = PDSpec(
        DiseaseProgressionPD(LogisticGrowth),
        "tumor_logistic",
        DiseaseProgressionPDParams(S0, kgrow, Smax, 0.0, 0.0),
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    tumor = res.observations[:tumor]

    # Logistic should approach Smax
    @test tumor[end] > S0
    @test tumor[end] < Smax

    # S-shaped growth: maximum growth rate at S = Smax/2
    @test tumor[1] < Smax / 2  # Start below inflection
end

@testset "Disease progression: Linear growth" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    S0 = 10.0
    alpha = 0.5  # Linear growth rate

    pd = PDSpec(
        DiseaseProgressionPD(LinearGrowth),
        "tumor_linear",
        DiseaseProgressionPDParams(S0, 0.0, 100.0, alpha, 0.0),  # kdrug=0
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    # Linear growth: S(t) = S0 + alpha * t
    for (i, t) in enumerate(res.t)
        expected_S = S0 + alpha * t
        @test isapprox(res.observations[:tumor][i], expected_S; rtol=1e-4)
    end
end

@testset "Disease progression: Drug effect reduces tumor" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(2.0, 50.0),  # Slower clearance = longer exposure
        [DoseEvent(0.0, 500.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    S0 = 50.0
    kgrow = 0.02
    Smax = 100.0
    kdrug = 0.01  # Drug effect

    # Without drug
    pd_no_drug = PDSpec(
        DiseaseProgressionPD(ExponentialGrowth),
        "tumor_no_drug",
        DiseaseProgressionPDParams(S0, kgrow, Smax, 0.0, 0.0),
        :conc,
        :tumor,
    )

    # With drug
    pd_with_drug = PDSpec(
        DiseaseProgressionPD(ExponentialGrowth),
        "tumor_with_drug",
        DiseaseProgressionPDParams(S0, kgrow, Smax, 0.0, kdrug),
        :conc,
        :tumor,
    )

    res_no_drug = simulate_pkpd_coupled(pk, pd_no_drug, grid, solver)
    res_with_drug = simulate_pkpd_coupled(pk, pd_with_drug, grid, solver)

    tumor_no_drug = res_no_drug.observations[:tumor]
    tumor_with_drug = res_with_drug.observations[:tumor]

    # Drug should reduce tumor growth
    @test tumor_with_drug[end] < tumor_no_drug[end]

    # Both should start at S0
    @test tumor_no_drug[1] ≈ S0
    @test tumor_with_drug[1] ≈ S0
end

@testset "Disease progression: Asymptotic growth" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 100.0, collect(0.0:2.0:100.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    S0 = 10.0
    kgrow = 0.05
    Smax = 80.0

    pd = PDSpec(
        DiseaseProgressionPD(AsymptoticGrowth),
        "tumor_asymptotic",
        DiseaseProgressionPDParams(S0, kgrow, Smax, 0.0, 0.0),
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    tumor = res.observations[:tumor]

    # Asymptotic growth: approaches Smax exponentially
    @test tumor[end] > S0
    @test tumor[end] < Smax

    # Growth rate decreases as S approaches Smax
    @test tumor[end] > 0.8 * Smax  # Should be close to Smax
end

@testset "Disease progression: metadata contains growth model" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        DiseaseProgressionPD(GompertzGrowth),
        "tumor_gompertz",
        DiseaseProgressionPDParams(10.0, 0.05, 100.0, 0.0, 0.0),
        :conc,
        :tumor,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test res.metadata["pd_kind"] == "DiseaseProgressionPD"
    @test res.metadata["growth_model"] == "GompertzGrowth"
end

# =============================================================================
# Combination Effect Models Tests
# Bliss Independence, Competitive Inhibition, Drug Interaction
# =============================================================================

@testset "Bliss Independence: basic evaluation" begin
    pd = PDSpec(
        BlissIndependence(),
        "bliss_combo",
        BlissIndependenceParams(
            0.0,      # E0
            1.0,      # Emax_A
            1.0,      # EC50_A
            1.0,      # gamma_A
            1.0,      # Emax_B
            1.0,      # EC50_B
            1.0,      # gamma_B
            :conc_A,  # input_A
            :conc_B,  # input_B
        ),
        :conc_A,
        :effect,
    )

    # Validate parameters
    validate(pd)

    # Test with no drug (both concentrations = 0)
    obs_zero = Dict(:conc_A => [0.0], :conc_B => [0.0])
    effect_zero = evaluate(pd, obs_zero)
    @test effect_zero[1] ≈ 0.0  # E0 = 0, no drug effect

    # Test with only drug A at EC50 (E_A = 0.5 * Emax_A = 0.5)
    obs_A_only = Dict(:conc_A => [1.0], :conc_B => [0.0])
    effect_A = evaluate(pd, obs_A_only)
    @test effect_A[1] ≈ 0.5  # E0 + E_A = 0 + 0.5

    # Test with only drug B at EC50 (E_B = 0.5)
    obs_B_only = Dict(:conc_A => [0.0], :conc_B => [1.0])
    effect_B = evaluate(pd, obs_B_only)
    @test effect_B[1] ≈ 0.5

    # Test Bliss Independence: E = E_A + E_B - E_A * E_B
    # With E_A = E_B = 0.5: E = 0.5 + 0.5 - 0.5*0.5 = 0.75
    obs_both = Dict(:conc_A => [1.0], :conc_B => [1.0])
    effect_both = evaluate(pd, obs_both)
    @test effect_both[1] ≈ 0.75
end

@testset "Bliss Independence: high concentrations approach maximum" begin
    pd = PDSpec(
        BlissIndependence(),
        "bliss_high",
        BlissIndependenceParams(
            0.0,      # E0
            1.0,      # Emax_A (fractional)
            1.0,      # EC50_A
            1.0,      # gamma_A
            1.0,      # Emax_B (fractional)
            1.0,      # EC50_B
            1.0,      # gamma_B
            :conc_A,
            :conc_B,
        ),
        :conc_A,
        :effect,
    )

    # Very high concentrations: both effects approach Emax (1.0)
    # Bliss: 1.0 + 1.0 - 1.0*1.0 = 1.0
    obs_high = Dict(:conc_A => [1000.0], :conc_B => [1000.0])
    effect_high = evaluate(pd, obs_high)
    @test effect_high[1] ≈ 1.0 atol=0.01
end

@testset "Bliss Independence: vector inputs" begin
    pd = PDSpec(
        BlissIndependence(),
        "bliss_vector",
        BlissIndependenceParams(5.0, 10.0, 1.0, 1.0, 10.0, 1.0, 1.0, :A, :B),
        :A,
        :effect,
    )

    # Test with time series
    obs = Dict(
        :A => [0.0, 0.5, 1.0, 2.0],
        :B => [0.0, 0.5, 1.0, 2.0],
    )

    effect = evaluate(pd, obs)
    @test length(effect) == 4

    # Effect should increase with concentration
    @test effect[1] < effect[2] < effect[3] < effect[4]

    # First point (no drug) should be E0
    @test effect[1] ≈ 5.0
end

@testset "Competitive Inhibition: basic evaluation" begin
    pd = PDSpec(
        CompetitiveInhibition(),
        "comp_inhib",
        CompetitiveInhibitionParams(
            0.0,        # E0
            100.0,      # Emax
            1.0,        # EC50
            1.0,        # gamma
            1.0,        # Ki
            :drug,
            :inhibitor,
        ),
        :drug,
        :effect,
    )

    validate(pd)

    # No drug, no inhibitor
    obs_zero = Dict(:drug => [0.0], :inhibitor => [0.0])
    effect_zero = evaluate(pd, obs_zero)
    @test effect_zero[1] ≈ 0.0  # E0

    # Drug at EC50, no inhibitor: E = E0 + Emax * 0.5 = 50
    obs_drug_only = Dict(:drug => [1.0], :inhibitor => [0.0])
    effect_drug = evaluate(pd, obs_drug_only)
    @test effect_drug[1] ≈ 50.0

    # Drug at EC50, inhibitor at Ki: EC50_apparent = EC50 * (1 + 1/1) = 2
    # Effect = Emax * C^gamma / (EC50_apparent^gamma + C^gamma)
    #        = 100 * 1 / (2 + 1) = 100/3 ≈ 33.33
    obs_with_inhib = Dict(:drug => [1.0], :inhibitor => [1.0])
    effect_inhib = evaluate(pd, obs_with_inhib)
    @test effect_inhib[1] ≈ 100.0 / 3.0 atol=0.1

    # Inhibitor reduces effect (competitive)
    @test effect_inhib[1] < effect_drug[1]
end

@testset "Competitive Inhibition: increasing inhibitor shifts EC50" begin
    pd = PDSpec(
        CompetitiveInhibition(),
        "comp_inhib_shift",
        CompetitiveInhibitionParams(0.0, 100.0, 1.0, 1.0, 0.5, :drug, :inhibitor),
        :drug,
        :effect,
    )

    # High drug concentration overcomes inhibition
    obs_high_drug = Dict(:drug => [100.0], :inhibitor => [1.0])
    effect_high = evaluate(pd, obs_high_drug)

    # Should still approach Emax at very high drug concentrations
    @test effect_high[1] > 90.0
end

@testset "Competitive Inhibition: vector inputs" begin
    pd = PDSpec(
        CompetitiveInhibition(),
        "comp_inhib_vec",
        CompetitiveInhibitionParams(0.0, 100.0, 1.0, 1.0, 1.0, :drug, :inhibitor),
        :drug,
        :effect,
    )

    # Varying drug, constant inhibitor
    obs = Dict(
        :drug => [0.0, 0.5, 1.0, 2.0, 5.0],
        :inhibitor => [1.0, 1.0, 1.0, 1.0, 1.0],
    )

    effect = evaluate(pd, obs)
    @test length(effect) == 5

    # Effect should increase with drug concentration (even with inhibitor)
    @test effect[1] < effect[2] < effect[3] < effect[4] < effect[5]
end

@testset "Drug Interaction (Greco): additive case (psi=0)" begin
    pd = PDSpec(
        DrugInteraction(),
        "greco_additive",
        DrugInteractionParams(
            0.0,      # E0
            100.0,    # Emax
            1.0,      # EC50_A
            1.0,      # EC50_B
            0.0,      # psi = 0 (additive)
            :A,
            :B,
        ),
        :A,
        :effect,
    )

    validate(pd)

    # No drugs
    obs_zero = Dict(:A => [0.0], :B => [0.0])
    effect_zero = evaluate(pd, obs_zero)
    @test effect_zero[1] ≈ 0.0

    # Drug A only at EC50: a=1, b=0, interaction=0
    # E = Emax * (1 + 0 + 0) / (1 + 1 + 0 + 0) = 100 * 1/2 = 50
    obs_A = Dict(:A => [1.0], :B => [0.0])
    effect_A = evaluate(pd, obs_A)
    @test effect_A[1] ≈ 50.0

    # Both drugs at EC50: a=1, b=1, interaction=0
    # E = Emax * (1 + 1 + 0) / (1 + 1 + 1 + 0) = 100 * 2/3 ≈ 66.67
    obs_both = Dict(:A => [1.0], :B => [1.0])
    effect_both = evaluate(pd, obs_both)
    @test effect_both[1] ≈ 100 * 2 / 3 atol=0.1
end

@testset "Drug Interaction (Greco): synergy case (psi>0)" begin
    pd = PDSpec(
        DrugInteraction(),
        "greco_synergy",
        DrugInteractionParams(
            0.0,      # E0
            100.0,    # Emax
            1.0,      # EC50_A
            1.0,      # EC50_B
            2.0,      # psi > 0 (synergy)
            :A,
            :B,
        ),
        :A,
        :effect,
    )

    # Both drugs at EC50: a=1, b=1, interaction=2*1*1=2
    # E = Emax * (1 + 1 + 2) / (1 + 1 + 1 + 2) = 100 * 4/5 = 80
    obs_both = Dict(:A => [1.0], :B => [1.0])
    effect_synergy = evaluate(pd, obs_both)
    @test effect_synergy[1] ≈ 80.0

    # Compare to additive: synergy should give higher effect
    pd_additive = PDSpec(
        DrugInteraction(),
        "greco_add",
        DrugInteractionParams(0.0, 100.0, 1.0, 1.0, 0.0, :A, :B),
        :A,
        :effect,
    )
    effect_additive = evaluate(pd_additive, obs_both)

    @test effect_synergy[1] > effect_additive[1]
end

@testset "Drug Interaction (Greco): antagonism case (psi<0)" begin
    pd = PDSpec(
        DrugInteraction(),
        "greco_antagonism",
        DrugInteractionParams(
            0.0,      # E0
            100.0,    # Emax
            1.0,      # EC50_A
            1.0,      # EC50_B
            -0.5,     # psi < 0 (antagonism)
            :A,
            :B,
        ),
        :A,
        :effect,
    )

    # Both drugs at EC50: a=1, b=1, interaction=-0.5*1*1=-0.5
    # E = Emax * (1 + 1 - 0.5) / (1 + 1 + 1 - 0.5) = 100 * 1.5/2.5 = 60
    obs_both = Dict(:A => [1.0], :B => [1.0])
    effect_antag = evaluate(pd, obs_both)
    @test effect_antag[1] ≈ 60.0

    # Compare to additive: antagonism should give lower effect
    pd_additive = PDSpec(
        DrugInteraction(),
        "greco_add",
        DrugInteractionParams(0.0, 100.0, 1.0, 1.0, 0.0, :A, :B),
        :A,
        :effect,
    )
    effect_additive = evaluate(pd_additive, obs_both)

    @test effect_antag[1] < effect_additive[1]
end

@testset "Drug Interaction (Greco): vector inputs" begin
    pd = PDSpec(
        DrugInteraction(),
        "greco_vec",
        DrugInteractionParams(0.0, 100.0, 1.0, 1.0, 1.0, :A, :B),
        :A,
        :effect,
    )

    # Varying both drugs together
    obs = Dict(
        :A => [0.0, 0.5, 1.0, 2.0],
        :B => [0.0, 0.5, 1.0, 2.0],
    )

    effect = evaluate(pd, obs)
    @test length(effect) == 4

    # Effect should increase with concentration
    @test effect[1] < effect[2] < effect[3] < effect[4]

    # First point should be E0
    @test effect[1] ≈ 0.0
end

@testset "Combination models: validation errors" begin
    # Bliss: negative EC50
    @test_throws Exception validate(PDSpec(
        BlissIndependence(),
        "bliss_bad",
        BlissIndependenceParams(0.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, :A, :B),
        :A,
        :effect,
    ))

    # Bliss: negative gamma
    @test_throws Exception validate(PDSpec(
        BlissIndependence(),
        "bliss_bad",
        BlissIndependenceParams(0.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, :A, :B),
        :A,
        :effect,
    ))

    # Competitive: negative Ki
    @test_throws Exception validate(PDSpec(
        CompetitiveInhibition(),
        "comp_bad",
        CompetitiveInhibitionParams(0.0, 100.0, 1.0, 1.0, -1.0, :drug, :inhibitor),
        :drug,
        :effect,
    ))

    # Drug Interaction: negative EC50
    @test_throws Exception validate(PDSpec(
        DrugInteraction(),
        "greco_bad",
        DrugInteractionParams(0.0, 100.0, -1.0, 1.0, 0.0, :A, :B),
        :A,
        :effect,
    ))
end

@testset "Combination models: missing input data errors" begin
    pd_bliss = PDSpec(
        BlissIndependence(),
        "bliss",
        BlissIndependenceParams(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, :A, :B),
        :A,
        :effect,
    )

    # Missing input_B
    @test_throws Exception evaluate(pd_bliss, Dict(:A => [1.0]))

    pd_comp = PDSpec(
        CompetitiveInhibition(),
        "comp",
        CompetitiveInhibitionParams(0.0, 100.0, 1.0, 1.0, 1.0, :drug, :inhibitor),
        :drug,
        :effect,
    )

    # Missing inhibitor
    @test_throws Exception evaluate(pd_comp, Dict(:drug => [1.0]))

    pd_greco = PDSpec(
        DrugInteraction(),
        "greco",
        DrugInteractionParams(0.0, 100.0, 1.0, 1.0, 0.0, :A, :B),
        :A,
        :effect,
    )

    # Missing input_B
    @test_throws Exception evaluate(pd_greco, Dict(:A => [1.0]))
end

# =============================================================================
# Tolerance Models Tests
# Counter-Regulation and Receptor Regulation
# =============================================================================

@testset "Tolerance Counter-Regulation: no tolerance with alpha=0" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    E0 = 10.0
    Emax = 50.0
    EC50 = 1.0
    gamma = 1.0

    pd = PDSpec(
        ToleranceCounterRegulation(),
        "tolerance_no_effect",
        ToleranceCounterRegulationParams(E0, Emax, EC50, gamma, 0.1, 0.1, 0.0),  # alpha=0
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test haskey(res.observations, :E_drug)
    @test haskey(res.states, :M)

    # With alpha=0, effect should equal E0 + E_drug (no tolerance)
    for i in eachindex(res.t)
        expected_effect = E0 + res.observations[:E_drug][i]
        @test isapprox(res.observations[:effect][i], expected_effect; rtol=1e-6)
    end
end

@testset "Tolerance Counter-Regulation: tolerance develops over time" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(2.0, 50.0),  # Slower clearance
        [DoseEvent(0.0, 300.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    E0 = 10.0
    Emax = 50.0
    EC50 = 1.0
    gamma = 1.0
    kin_mod = 0.2  # Moderator accumulation rate
    kout_mod = 0.05  # Slow moderator elimination
    alpha = 1.0  # Strong feedback

    pd = PDSpec(
        ToleranceCounterRegulation(),
        "tolerance_develops",
        ToleranceCounterRegulationParams(E0, Emax, EC50, gamma, kin_mod, kout_mod, alpha),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    effect = res.observations[:effect]
    E_drug = res.observations[:E_drug]
    moderator = res.states[:M]

    # Moderator should start at 0
    @test isapprox(moderator[1], 0.0; atol=1e-6)

    # Moderator should increase over time with drug exposure
    @test moderator[end] > moderator[1]

    # Check tolerance at a later time point (after moderator has accumulated)
    # At t=10 (index 11), there should be significant drug effect and some tolerance
    mid_idx = 11  # t=10

    # Effect with tolerance should be less than effect without tolerance
    # E_net = E0 + E_drug - alpha * M, so if M > 0, E_net < E0 + E_drug
    effect_without_tolerance = E0 + E_drug[mid_idx]
    @test effect[mid_idx] < effect_without_tolerance
end

@testset "Tolerance Counter-Regulation: repeated dosing builds tolerance" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0), DoseEvent(48.0, 100.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        ToleranceCounterRegulation(),
        "tolerance_repeated",
        ToleranceCounterRegulationParams(10.0, 50.0, 1.0, 1.0, 0.2, 0.05, 1.0),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    effect = res.observations[:effect]
    conc = res.observations[:conc]

    # Find peak effects after each dose
    # Dose 1: around t=0-1
    peak1_idx = argmax(effect[1:10])

    # Dose 2: around t=24-25
    peak2_idx = 24 + argmax(effect[25:35])

    # Dose 3: around t=48-49
    peak3_idx = 48 + argmax(effect[49:60])

    # Tolerance should build: peak effect should decrease with each dose
    # (assuming similar concentrations at peaks)
    @test effect[peak1_idx] > effect[peak2_idx]  # Tolerance develops
end

@testset "Receptor Regulation Down: receptor decreases with drug exposure" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(2.0, 50.0),
        [DoseEvent(0.0, 300.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        ReceptorRegulation(),
        "receptor_down",
        ReceptorRegulationParams(
            10.0,   # E0
            50.0,   # Emax
            1.0,    # EC50
            1.0,    # gamma
            1.0,    # R_baseline (normalized)
            0.05,   # kreg (slow return to baseline)
            2.0,    # Rmax (for up-regulation)
            0.2,    # kchange (rate of down-regulation)
            :down,  # direction
        ),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    @test haskey(res.observations, :receptor_density)
    @test haskey(res.states, :R)

    receptor = res.states[:R]

    # Receptor should start at baseline (1.0)
    @test isapprox(receptor[1], 1.0; atol=1e-6)

    # With down-regulation, receptor should decrease below baseline
    @test minimum(receptor) < 1.0

    # Effect should be reduced due to lower receptor density
    effect = res.observations[:effect]
    E0 = 10.0

    # At some point, effect should be less than if receptors were at baseline
    @test any(e < E0 + 50.0 for e in effect[2:end])
end

@testset "Receptor Regulation Up: receptor increases with drug exposure" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(2.0, 50.0),
        [DoseEvent(0.0, 300.0)],
    )

    grid = SimGrid(0.0, 72.0, collect(0.0:1.0:72.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        ReceptorRegulation(),
        "receptor_up",
        ReceptorRegulationParams(
            10.0,   # E0
            50.0,   # Emax
            1.0,    # EC50
            1.0,    # gamma
            1.0,    # R_baseline
            0.05,   # kreg
            2.0,    # Rmax
            0.2,    # kchange (rate of up-regulation)
            :up,    # direction
        ),
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    receptor = res.states[:R]

    # Receptor should start at baseline (1.0)
    @test isapprox(receptor[1], 1.0; atol=1e-6)

    # With up-regulation, receptor should increase above baseline
    @test maximum(receptor) > 1.0

    # Receptor should not exceed Rmax
    @test all(r <= 2.0 + 1e-6 for r in receptor)
end

@testset "Receptor Regulation: no change with kchange=0" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 200.0)],
    )

    grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    pd = PDSpec(
        ReceptorRegulation(),
        "receptor_stable",
        ReceptorRegulationParams(10.0, 50.0, 1.0, 1.0, 1.0, 0.1, 2.0, 0.0, :down),  # kchange=0
        :conc,
        :effect,
    )

    res = simulate_pkpd_coupled(pk, pd, grid, solver)

    receptor = res.states[:R]

    # With kchange=0, receptor should stay at baseline
    for r in receptor
        @test isapprox(r, 1.0; rtol=1e-6)
    end
end

@testset "Tolerance models: validation errors" begin
    # Counter-regulation: negative EC50
    @test_throws Exception validate(PDSpec(
        ToleranceCounterRegulation(),
        "bad",
        ToleranceCounterRegulationParams(10.0, 50.0, -1.0, 1.0, 0.1, 0.1, 1.0),
        :conc,
        :effect,
    ))

    # Counter-regulation: negative kin_mod
    @test_throws Exception validate(PDSpec(
        ToleranceCounterRegulation(),
        "bad",
        ToleranceCounterRegulationParams(10.0, 50.0, 1.0, 1.0, -0.1, 0.1, 1.0),
        :conc,
        :effect,
    ))

    # Counter-regulation: negative alpha
    @test_throws Exception validate(PDSpec(
        ToleranceCounterRegulation(),
        "bad",
        ToleranceCounterRegulationParams(10.0, 50.0, 1.0, 1.0, 0.1, 0.1, -1.0),
        :conc,
        :effect,
    ))

    # Receptor regulation: invalid direction
    @test_throws Exception validate(PDSpec(
        ReceptorRegulation(),
        "bad",
        ReceptorRegulationParams(10.0, 50.0, 1.0, 1.0, 1.0, 0.1, 2.0, 0.1, :invalid),
        :conc,
        :effect,
    ))

    # Receptor regulation: negative R_baseline
    @test_throws Exception validate(PDSpec(
        ReceptorRegulation(),
        "bad",
        ReceptorRegulationParams(10.0, 50.0, 1.0, 1.0, -1.0, 0.1, 2.0, 0.1, :down),
        :conc,
        :effect,
    ))
end

@testset "Tolerance models: metadata correct" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    # Counter-regulation
    pd_tol = PDSpec(
        ToleranceCounterRegulation(),
        "tolerance",
        ToleranceCounterRegulationParams(10.0, 50.0, 1.0, 1.0, 0.1, 0.1, 1.0),
        :conc,
        :effect,
    )

    res_tol = simulate_pkpd_coupled(pk, pd_tol, grid, solver)
    @test res_tol.metadata["pd_kind"] == "ToleranceCounterRegulation"

    # Receptor regulation
    pd_rec = PDSpec(
        ReceptorRegulation(),
        "receptor",
        ReceptorRegulationParams(10.0, 50.0, 1.0, 1.0, 1.0, 0.1, 2.0, 0.1, :down),
        :conc,
        :effect,
    )

    res_rec = simulate_pkpd_coupled(pk, pd_rec, grid, solver)
    @test res_rec.metadata["pd_kind"] == "ReceptorRegulation"
    @test res_rec.metadata["direction"] == "down"
end
