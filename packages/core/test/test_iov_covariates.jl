# IOV and covariate tests

@testset "IOV changes parameters across occasions deterministically" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "iov_pop_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    iiv = nothing
    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(1234), OccasionDefinition(:dose_times)
    )

    pop = PopulationSpec(base, iiv, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res1 = simulate_population(pop, grid, solver)
    res2 = simulate_population(pop, grid, solver)

    # Deterministic replay
    c1 = res1.individuals[1].observations[:conc]
    c2 = res2.individuals[1].observations[:conc]
    for i in eachindex(c1)
        @test isapprox(c1[i], c2[i]; rtol=1e-12, atol=1e-12)
    end

    # Behavioral: after second dose, concentration trajectory should differ from a non-IOV run
    pop_no_iov = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_no_iov, grid, solver)

    c_base = base_res.individuals[1].observations[:conc]
    # pick a time after 12.0, for example index where t == 13
    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(c1[idx13] - c_base[idx13]) > 0.0
end

@testset "Coupled PKPD with IOV is deterministic and alters post-dose dynamics" begin
    pk = ModelSpec(
        OneCompIVBolus(),
        "pkpd_iov",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0), DoseEvent(12.0, 100.0)],
    )

    Kin = 10.0
    Kout = 0.5
    Rss = Kin / Kout

    pd = PDSpec(
        IndirectResponseTurnover(),
        "pd_turnover",
        IndirectResponseTurnoverParams(Kin, Kout, Rss, 0.8, 0.5),
        :conc,
        :response,
    )

    iov = IOVSpec(
        LogNormalIIV(), Dict(:CL => 0.3), UInt64(2222), OccasionDefinition(:dose_times)
    )
    pop = PopulationSpec(pk, nothing, iov, nothing, IndividualCovariates[])

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    r1 = simulate_population(pop, grid, solver; pd_spec=pd)
    r2 = simulate_population(pop, grid, solver; pd_spec=pd)

    # Deterministic replay
    y1 = r1.individuals[1].observations[:response]
    y2 = r2.individuals[1].observations[:response]

    for i in eachindex(y1)
        @test isapprox(y1[i], y2[i]; rtol=1e-12, atol=1e-12)
    end

    # Compare against no-IOV behavior
    pop_no_iov = PopulationSpec(pk, nothing, nothing, nothing, IndividualCovariates[])
    base = simulate_population(pop_no_iov, grid, solver; pd_spec=pd)

    idx13 = findfirst(==(13.0), grid.saveat)
    @test idx13 !== nothing
    @test abs(y1[idx13] - base.individuals[1].observations[:response][idx13]) > 0.0
end

@testset "Time-varying covariate changes concentration deterministically via segmentation" begin
    base = ModelSpec(
        OneCompIVBolus(),
        "tv_cov_iv",
        OneCompIVBolusParams(5.0, 50.0),
        [DoseEvent(0.0, 100.0)],
    )

    # CL increases at t=10, so concentrations after 10 should drop faster
    cm = CovariateModel(
        "cl_tv", [CovariateEffect(LinearCovariate(), :CL, :CLMULT, 1.0, 1.0)]
    )

    tv = TimeVaryingCovariates(
        Dict(:CLMULT => TimeCovariateSeries(StepTimeCovariate(), [0.0, 10.0], [1.0, 2.0]))
    )

    covs = [IndividualCovariates(Dict{Symbol,Float64}(), tv)]
    pop = PopulationSpec(base, nothing, nothing, cm, covs)

    grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
    solver = SolverSpec(:Tsit5, 1e-10, 1e-12, 10^7)

    res = simulate_population(pop, grid, solver)

    c = res.individuals[1].observations[:conc]
    idx9 = findfirst(==(9.0), grid.saveat)
    idx11 = findfirst(==(11.0), grid.saveat)

    @test idx9 !== nothing
    @test idx11 !== nothing

    # sanity: conc at 11 should be lower than it would be with constant CL
    pop_const = PopulationSpec(base, nothing, nothing, nothing, IndividualCovariates[])
    base_res = simulate_population(pop_const, grid, solver)
    c_base = base_res.individuals[1].observations[:conc]

    @test abs(c[idx11] - c_base[idx11]) > 0.0
end
