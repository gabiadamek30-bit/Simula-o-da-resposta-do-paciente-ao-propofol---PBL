# Test suite for Covariate-on-IIV and IOV integration in FOCE-I

using Test
using NeoPKPD
using LinearAlgebra

@testset "Covariate-IIV and IOV Integration" begin

    @testset "CovariateOnIIV Type" begin
        # Test basic construction
        cov_effect = CovariateOnIIV(:eta_CL, :WT, 70.0)
        @test cov_effect.eta_name == :eta_CL
        @test cov_effect.covariate_name == :WT
        @test cov_effect.reference_value == 70.0
        @test cov_effect.effect_type == :exponential

        # Test linear effect type
        cov_linear = CovariateOnIIV(:eta_V, :AGE, 40.0; effect_type=:linear)
        @test cov_linear.effect_type == :linear
    end

    @testset "EstimationIOVSpec Type" begin
        iov_spec = EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2])
        @test iov_spec.eta_name == :eta_CL
        @test length(iov_spec.occasion_names) == 2
        @test iov_spec.omega_iov == 0.04  # Default

        # Custom omega_iov
        iov_spec2 = EstimationIOVSpec(:eta_V, [:V1, :V2, :V3]; omega_iov=0.09)
        @test iov_spec2.omega_iov == 0.09
    end

    @testset "OccasionData Type" begin
        occ_data = OccasionData("S001", [1, 1, 1, 2, 2, 2])
        @test occ_data.subject_id == "S001"
        @test length(occ_data.occasion_indices) == 6
    end

    @testset "compute_covariate_adjusted_omega" begin
        omega_base = diagm([0.09, 0.04])
        covariate_effects = [CovariateOnIIV(:eta_CL, :WT, 70.0; effect_type=:exponential)]

        # Subject with reference weight - no adjustment
        covariates_ref = Dict{Symbol,Float64}(:WT => 70.0)
        theta_cov = [0.0]  # No effect
        omega_names = [:eta_CL, :eta_V]

        omega_adj = compute_covariate_adjusted_omega(
            omega_base, covariate_effects, covariates_ref, theta_cov, omega_names
        )
        @test omega_adj ≈ omega_base

        # Subject with different weight - should adjust
        covariates_heavy = Dict{Symbol,Float64}(:WT => 90.0)
        theta_cov2 = [0.5]  # Positive effect
        omega_adj2 = compute_covariate_adjusted_omega(
            omega_base, covariate_effects, covariates_heavy, theta_cov2, omega_names
        )
        # Omega for CL should be higher for heavier subjects
        @test omega_adj2[1,1] > omega_base[1,1]
        @test omega_adj2[2,2] ≈ omega_base[2,2]  # V not affected
    end

    @testset "compute_iov_eta" begin
        eta_iiv = [0.1, 0.05]  # IIV for CL and V
        iov_specs = [EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2]; omega_iov=0.04)]
        eta_iov = [0.02 -0.02; 0.0 0.0]  # IOV etas: 2 params × 2 occasions
        omega_names = [:eta_CL, :eta_V]

        # Occasion 1
        eta_total_1 = compute_iov_eta(eta_iiv, iov_specs, 1, eta_iov, omega_names)
        @test eta_total_1[1] ≈ eta_iiv[1] + eta_iov[1, 1]
        @test eta_total_1[2] ≈ eta_iiv[2]  # V has no IOV

        # Occasion 2
        eta_total_2 = compute_iov_eta(eta_iiv, iov_specs, 2, eta_iov, omega_names)
        @test eta_total_2[1] ≈ eta_iiv[1] + eta_iov[1, 2]
    end

    @testset "get_n_iov_params" begin
        iov_specs_empty = EstimationIOVSpec[]
        @test get_n_iov_params(iov_specs_empty) == 0

        iov_specs_one = [EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2])]
        @test get_n_iov_params(iov_specs_one) == 2

        iov_specs_multi = [
            EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2]),
            EstimationIOVSpec(:eta_V, [:OCC1, :OCC2, :OCC3])
        ]
        @test get_n_iov_params(iov_specs_multi) == 5
    end

    @testset "build_iov_omega" begin
        iov_specs = [
            EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2]; omega_iov=0.04),
            EstimationIOVSpec(:eta_V, [:OCC1, :OCC2]; omega_iov=0.01)
        ]

        omega_iov = build_iov_omega(iov_specs)
        @test size(omega_iov) == (4, 4)
        @test omega_iov[1,1] ≈ 0.04
        @test omega_iov[2,2] ≈ 0.04
        @test omega_iov[3,3] ≈ 0.01
        @test omega_iov[4,4] ≈ 0.01
        @test omega_iov[1,2] ≈ 0.0  # Off-diagonal should be zero
    end

    @testset "extract_occasion_indices" begin
        # With OCC key
        covariates_occ = Dict{Symbol,Any}(:OCC => [1, 1, 2, 2])
        indices = extract_occasion_indices(covariates_occ, 4)
        @test indices == [1, 1, 2, 2]

        # With OCCASION key
        covariates_occasion = Dict{Symbol,Any}(:OCCASION => [1, 2, 3])
        indices2 = extract_occasion_indices(covariates_occasion, 3)
        @test indices2 == [1, 2, 3]

        # No occasion data - defaults to all 1s
        covariates_none = Dict{Symbol,Any}(:WT => 70.0)
        indices3 = extract_occasion_indices(covariates_none, 5)
        @test indices3 == [1, 1, 1, 1, 1]
    end

    @testset "FOCE-I with Covariate-on-IIV" begin
        doses = [DoseEvent(0.0, 100.0)]

        # Create subjects with different weights
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses;
            covariates=Dict{Symbol,Any}(:WT => 70.0))
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses;
            covariates=Dict{Symbol,Any}(:WT => 85.0))
        subj3 = SubjectData("S3", [0.5, 1.0, 2.0, 4.0], [1.5, 1.3, 1.1, 0.7], doses;
            covariates=Dict{Symbol,Any}(:WT => 55.0))
        observed = ObservedData([subj1, subj2, subj3])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Specify covariate effect on IIV for CL
        wt_effect = CovariateOnIIV(:eta_CL, :WT, 70.0; effect_type=:exponential)

        config = EstimationConfig(
            FOCEIMethod(max_inner_iter=20, inner_tol=1e-4);
            theta_init=[10.0, 50.0],
            theta_lower=[1.0, 10.0],
            theta_upper=[100.0, 200.0],
            theta_names=[:CL, :V],
            omega_init=diagm([0.09, 0.04]),
            omega_names=[:eta_CL, :eta_V],
            sigma_init=sigma_spec,
            max_iter=15,
            tol=1e-2,
            compute_se=false,
            verbose=false,
            covariate_on_iiv=[wt_effect]
        )

        result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test isfinite(result.ofv)
        @test result.theta[1] > 0  # CL > 0
        @test result.theta[2] > 0  # V > 0
        @test any(contains.(result.messages, "Covariate effects on IIV"))
    end

    @testset "FOCE-I with IOV specification" begin
        doses = [DoseEvent(0.0, 100.0)]

        # Create subjects with occasion data
        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses;
            covariates=Dict{Symbol,Any}(:OCC => [1, 1, 1, 1]))
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses;
            covariates=Dict{Symbol,Any}(:OCC => [1, 1, 1, 1]))
        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        # Specify IOV on CL
        iov_cl = EstimationIOVSpec(:eta_CL, [:OCC1, :OCC2]; omega_iov=0.04)

        config = EstimationConfig(
            FOCEIMethod(max_inner_iter=20, inner_tol=1e-4);
            theta_init=[10.0, 50.0],
            theta_lower=[1.0, 10.0],
            theta_upper=[100.0, 200.0],
            theta_names=[:CL, :V],
            omega_init=diagm([0.09, 0.04]),
            omega_names=[:eta_CL, :eta_V],
            sigma_init=sigma_spec,
            max_iter=10,
            tol=1e-2,
            compute_se=false,
            verbose=false,
            iov_specs=[iov_cl]
        )

        result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test isfinite(result.ofv)
        @test any(contains.(result.messages, "IOV enabled"))
    end

    @testset "FOCE-I backward compatibility (no advanced features)" begin
        # Ensure basic estimation still works when no covariates or IOV specified
        doses = [DoseEvent(0.0, 100.0)]

        subj1 = SubjectData("S1", [0.5, 1.0, 2.0, 4.0], [1.8, 1.6, 1.3, 0.9], doses)
        subj2 = SubjectData("S2", [0.5, 1.0, 2.0, 4.0], [2.2, 1.9, 1.5, 1.0], doses)
        observed = ObservedData([subj1, subj2])

        model_params = OneCompIVBolusParams(10.0, 50.0)
        model_spec = ModelSpec(OneCompIVBolus(), "pk", model_params, doses)
        grid = SimGrid(0.0, 8.0, collect(0.0:0.5:8.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sigma_spec = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        config = EstimationConfig(
            FOCEIMethod();
            theta_init=[10.0, 50.0],
            theta_lower=[1.0, 10.0],
            theta_upper=[100.0, 200.0],
            omega_init=diagm([0.09, 0.04]),
            sigma_init=sigma_spec,
            max_iter=10,
            compute_se=false
        )

        result = NeoPKPD.estimate(observed, model_spec, config; grid=grid, solver=solver)

        @test result isa EstimationResult
        @test isfinite(result.ofv)
        @test !any(contains.(result.messages, "Covariate"))
        @test !any(contains.(result.messages, "IOV"))
    end

end
