# Sigmoid Emax PD model tests

@testset "SigmoidEmax basic evaluation" begin
    pd = PDSpec(
        SigmoidEmax(),
        "semax",
        SigmoidEmaxParams(10.0, 50.0, 1.0, 2.0),  # E0=10, Emax=50, EC50=1, gamma=2
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    effects = evaluate(pd, concentrations)

    for (i, C) in enumerate(concentrations)
        e_ref = sigmoid_emax_ref(C, 10.0, 50.0, 1.0, 2.0)
        @test isapprox(effects[i], e_ref; rtol=1e-10)
    end
end

@testset "SigmoidEmax: gamma=1 matches DirectEmax" begin
    pd_sigmoid = PDSpec(
        SigmoidEmax(),
        "semax_g1",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 1.0),  # gamma=1
        :conc,
        :effect,
    )

    pd_direct = PDSpec(
        DirectEmax(),
        "demax",
        DirectEmaxParams(10.0, 40.0, 0.8),
        :conc,
        :effect,
    )

    concentrations = [0.0, 0.1, 0.5, 1.0, 2.0, 5.0]

    effects_sig = evaluate(pd_sigmoid, concentrations)
    effects_dir = evaluate(pd_direct, concentrations)

    for i in eachindex(concentrations)
        @test isapprox(effects_sig[i], effects_dir[i]; rtol=1e-10)
    end
end

@testset "SigmoidEmax: higher gamma gives steeper response" begin
    pd_low = PDSpec(
        SigmoidEmax(),
        "semax_low",
        SigmoidEmaxParams(0.0, 100.0, 5.0, 1.0),  # gamma=1
        :conc,
        :effect,
    )

    pd_high = PDSpec(
        SigmoidEmax(),
        "semax_high",
        SigmoidEmaxParams(0.0, 100.0, 5.0, 4.0),  # gamma=4
        :conc,
        :effect,
    )

    # At EC50, both should give ~50% of Emax
    effects_low = evaluate(pd_low, [5.0])
    effects_high = evaluate(pd_high, [5.0])
    @test isapprox(effects_low[1], 50.0; rtol=0.01)
    @test isapprox(effects_high[1], 50.0; rtol=0.01)

    # Below EC50, high gamma gives lower effect
    effects_low_sub = evaluate(pd_low, [2.0])
    effects_high_sub = evaluate(pd_high, [2.0])
    @test effects_high_sub[1] < effects_low_sub[1]

    # Above EC50, high gamma gives higher effect
    effects_low_sup = evaluate(pd_low, [10.0])
    effects_high_sup = evaluate(pd_high, [10.0])
    @test effects_high_sup[1] > effects_low_sup[1]
end

@testset "SigmoidEmax PKPD coupling" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pk_semax",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    pd = PDSpec(
        SigmoidEmax(),
        "semax",
        SigmoidEmaxParams(10.0, 40.0, 0.8, 2.0),
        :conc,
        :effect,
    )

    res = simulate_pkpd(pk, pd, grid, solver)

    @test haskey(res.observations, :conc)
    @test haskey(res.observations, :effect)
    @test length(res.observations[:effect]) == length(res.t)

    # Verify effects match expected values
    for (i, t) in enumerate(res.t)
        C = res.observations[:conc][i]
        e_ref = sigmoid_emax_ref(C, 10.0, 40.0, 0.8, 2.0)
        @test isapprox(res.observations[:effect][i], e_ref; rtol=1e-10)
    end
end
