# One-compartment PK model tests

@testset "OneCompIVBolus analytic equivalence" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Dose event at nonzero time" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_delayed",
        OneCompIVBolusParams(3.0, 30.0),
        [DoseEvent(2.0, 120.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:0.25:12.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue  # skip discontinuity (left-continuous solver output)
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "Multiple bolus schedule" begin
    spec = ModelSpec(
        OneCompIVBolus(),
        "1c_iv_bolus_multi",
        OneCompIVBolusParams(4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(6.0, 50.0), DoseEvent(10.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        if any(d.time == t for d in spec.doses)
            continue
        end
        c_ref = analytic_onecomp_ivbolus_conc(t, spec.doses, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "OneCompOralFirstOrder analytic equivalence" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo",
        OneCompOralFirstOrderParams(1.2, 5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "OneCompOralFirstOrder multiple doses" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "1c_oral_fo_multi",
        OneCompOralFirstOrderParams(0.9, 4.0, 40.0),
        [DoseEvent(0.0, 80.0), DoseEvent(8.0, 50.0), DoseEvent(16.0, 50.0)],
    )

    grid = SimGrid(0.0, 24.0, collect(0.0:0.25:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate(spec, grid, solver)

    Ka = spec.params.Ka
    CL = spec.params.CL
    V = spec.params.V

    for (i, t) in enumerate(res.t)
        c_ref = analytic_onecomp_oral_first_order_conc(t, spec.doses, Ka, CL, V)
        @test isapprox(res.observations[:conc][i], c_ref; rtol=1e-8, atol=1e-10)
    end
end

@testset "State outputs are present and aligned" begin
    spec = ModelSpec(
        OneCompOralFirstOrder(),
        "state_presence",
        OneCompOralFirstOrderParams(1.0, 3.0, 30.0),
        [DoseEvent(0.0, 100.0)],
    )

    grid = SimGrid(0.0, 12.0, collect(0.0:1.0:12.0))
    solver = SolverSpec(:Tsit5, 1e-9, 1e-11, 10^7)

    res = simulate(spec, grid, solver)

    @test haskey(res.states, :A_gut)
    @test haskey(res.states, :A_central)
    @test haskey(res.observations, :conc)

    for v in values(res.states)
        @test length(v) == length(res.t)
    end
end
