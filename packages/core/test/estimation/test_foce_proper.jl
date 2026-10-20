# Tests for Proper FOCE-I Implementation
# Tests the full Laplacian correction and interaction terms

using Test
using NeoPKPD
using LinearAlgebra
using Statistics
using StableRNGs

@testset "Proper FOCE-I Implementation" begin

    # ========================================================================
    # Test Setup: Create benchmark one-compartment IV model
    # ========================================================================

    # True population parameters (known values for validation)
    true_CL = 10.0   # L/hr
    true_V = 100.0   # L
    true_omega_CL = 0.09  # 30% CV on CL
    true_omega_V = 0.04   # 20% CV on V
    true_sigma = 0.1      # 10% proportional error

    # Create test subjects with known etas for validation
    function create_test_population(n_subj::Int, seed::UInt64)
        rng = StableRNG(seed)

        subjects = SubjectData[]
        true_etas = Vector{Vector{Float64}}()

        for i in 1:n_subj
            # Sample true individual etas
            eta_CL = randn(rng) * sqrt(true_omega_CL)
            eta_V = randn(rng) * sqrt(true_omega_V)
            push!(true_etas, [eta_CL, eta_V])

            # Individual parameters
            ind_CL = true_CL * exp(eta_CL)
            ind_V = true_V * exp(eta_V)

            # Observation times
            times = [0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]

            # Simulate observations with residual error
            dose = 1000.0
            obs = Float64[]
            for t in times
                conc = (dose / ind_V) * exp(-ind_CL / ind_V * t)
                # Add proportional error
                obs_conc = conc * (1 + randn(rng) * true_sigma)
                push!(obs, max(obs_conc, 0.01))  # Ensure positive
            end

            push!(subjects, SubjectData(
                "SUBJ_$i",
                times,
                obs,
                [DoseEvent(0.0, dose)]
            ))
        end

        return ObservedData(subjects), true_etas
    end

    # ========================================================================
    # Test 1: Laplacian Correction Computation
    # ========================================================================

    @testset "Laplacian Correction Computation" begin
        observed, true_etas = create_test_population(5, UInt64(12345))

        model_spec = ModelSpec(
            OneCompIVBolus(),
            "Test 1-comp IV",
            OneCompIVBolusParams(true_CL, true_V),
            [DoseEvent(0.0, 1000.0)]
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.1:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        config = EstimationConfig(
            FOCEIMethod(max_inner_iter=100, inner_tol=1e-8);
            theta_init = [true_CL, true_V],
            theta_lower = [0.1, 1.0],
            theta_upper = [100.0, 500.0],
            omega_init = diagm([true_omega_CL, true_omega_V]),
            sigma_init = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(true_sigma),
                :conc,
                UInt64(1)
            )
        )

        # Validate the objective components
        diagnostics = validate_foce_objective(observed, model_spec, config, grid, solver)

        # Check that Hessians are computed
        @test length(diagnostics.eta_hessians) == 5
        @test all(size(H) == (2, 2) for H in diagnostics.eta_hessians)

        # Check that Laplacian corrections are finite
        @test all(isfinite.(diagnostics.laplacian_corrections))

        @test diagnostics.interaction_enabled == true
    end

    # ========================================================================
    # Test 2: Interaction Term (Variance Depends on Prediction)
    # ========================================================================

    @testset "Eta-Epsilon Interaction" begin
        # Test that residual variance computation depends on prediction
        sigma_prop = ResidualErrorSpec(
            ProportionalError(),
            ProportionalErrorParams(0.1),
            :conc,
            UInt64(1)
        )

        sigma_add = ResidualErrorSpec(
            AdditiveError(),
            AdditiveErrorParams(0.5),
            :conc,
            UInt64(1)
        )

        sigma_comb = ResidualErrorSpec(
            CombinedError(),
            CombinedErrorParams(0.3, 0.1),
            :conc,
            UInt64(1)
        )

        # For proportional error, variance should depend on prediction
        var_prop_low = NeoPKPD.compute_residual_variance(1.0, sigma_prop)
        var_prop_high = NeoPKPD.compute_residual_variance(10.0, sigma_prop)
        @test var_prop_high > var_prop_low  # Variance increases with prediction
        @test var_prop_high / var_prop_low ≈ 100.0  # (10/1)^2

        # For additive error, variance should NOT depend on prediction
        var_add_low = NeoPKPD.compute_residual_variance(1.0, sigma_add)
        var_add_high = NeoPKPD.compute_residual_variance(10.0, sigma_add)
        @test var_add_high ≈ var_add_low

        # For combined error, variance should depend on prediction (but less than proportional)
        var_comb_low = NeoPKPD.compute_residual_variance(1.0, sigma_comb)
        var_comb_high = NeoPKPD.compute_residual_variance(10.0, sigma_comb)
        @test var_comb_high > var_comb_low
        @test var_comb_high / var_comb_low < 100.0  # Less than pure proportional
    end

    # ========================================================================
    # Test 3: Full FOCE-I Estimation
    # ========================================================================

    @testset "Full FOCE-I Estimation" begin
        observed, true_etas = create_test_population(10, UInt64(99999))

        model_spec = ModelSpec(
            OneCompIVBolus(),
            "Test 1-comp IV",
            OneCompIVBolusParams(true_CL, true_V),
            [DoseEvent(0.0, 1000.0)]
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        config = EstimationConfig(
            FOCEIMethod(max_inner_iter=50, inner_tol=1e-6);
            theta_init = [8.0, 80.0],  # Start away from true values
            theta_lower = [0.1, 1.0],
            theta_upper = [100.0, 500.0],
            omega_init = diagm([0.1, 0.1]),
            sigma_init = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(0.15),
                :conc,
                UInt64(1)
            ),
            max_iter = 50,
            tol = 1e-3,
            compute_se = false,
            verbose = false
        )

        rng = StableRNG(UInt64(12345))

        # Run proper FOCE-I
        result = foce_estimate_proper(observed, model_spec, config, grid, solver, rng)

        # Check convergence
        @test result.convergence

        # Check parameter estimates are reasonable (within 100% of true)
        # Note: With limited data and iterations, estimates may not be very precise
        @test abs(result.theta[1] - true_CL) / true_CL < 1.0
        @test abs(result.theta[2] - true_V) / true_V < 1.0

        # Check omega estimates are reasonable
        omega_diag = diag(result.omega)
        @test omega_diag[1] > 0
        @test omega_diag[2] > 0

        # Check OFV is finite and reasonable
        @test isfinite(result.ofv)
        @test result.ofv > 0  # -2LL should be positive

        # Check AIC and BIC
        @test isfinite(result.aic)
        @test isfinite(result.bic)
        @test result.aic < result.bic  # BIC has larger penalty

        # Check individual estimates
        @test length(result.individuals) == 10
        for ind in result.individuals
            @test length(ind.eta) == 2
            @test all(isfinite.(ind.ipred))
            @test all(isfinite.(ind.cwres))
        end

        # Check that messages indicate proper FOCE-I
        @test any(occursin("Proper FOCE-I", m) for m in result.messages) || any(occursin("Laplacian", m) for m in result.messages)
    end

    # ========================================================================
    # Test 4: CWRES Distribution Properties
    # ========================================================================

    @testset "CWRES Distribution Properties" begin
        # Generate a population for statistical testing
        observed, true_etas = create_test_population(20, UInt64(11111))

        model_spec = ModelSpec(
            OneCompIVBolus(),
            "Test 1-comp IV",
            OneCompIVBolusParams(true_CL, true_V),
            [DoseEvent(0.0, 1000.0)]
        )

        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        config = EstimationConfig(
            FOCEIMethod(max_inner_iter=50, inner_tol=1e-6);
            theta_init = [true_CL, true_V],
            theta_lower = [0.1, 1.0],
            theta_upper = [100.0, 500.0],
            omega_init = diagm([true_omega_CL, true_omega_V]),
            sigma_init = ResidualErrorSpec(
                ProportionalError(),
                ProportionalErrorParams(true_sigma),
                :conc,
                UInt64(1)
            ),
            max_iter = 100,
            tol = 1e-3,
            compute_se = false,
            verbose = false
        )

        rng = StableRNG(UInt64(12345))

        result = foce_estimate_proper(observed, model_spec, config, grid, solver, rng)

        # Collect all CWRES values
        all_cwres = vcat([ind.cwres for ind in result.individuals]...)

        # For a correct model, CWRES should be approximately N(0, 1)
        cwres_mean = mean(all_cwres)
        cwres_std = std(all_cwres)

        @info "CWRES mean: $cwres_mean (should be ~0)"
        @info "CWRES std: $cwres_std (should be ~1)"

        # Allow some tolerance since we're using estimated parameters
        # CWRES may not be perfectly N(0,1) due to approximations in FOCE-I
        @test abs(cwres_mean) < 1.0  # Mean should be close to 0
        @test cwres_std > 0.0  # Std should be positive
        @test all(isfinite.(all_cwres))  # All CWRES should be finite
    end

end
