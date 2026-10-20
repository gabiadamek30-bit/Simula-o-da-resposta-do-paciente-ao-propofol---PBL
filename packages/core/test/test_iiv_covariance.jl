# Test suite for Population IIV Covariance
# Tests full covariance matrix support with correlations

using Test
using NeoPKPD
using LinearAlgebra
using StableRNGs
using Statistics

@testset "IIV Covariance" begin

    @testset "OmegaMatrix Construction" begin
        @testset "Diagonal OmegaMatrix from Dict" begin
            omegas = Dict(:CL => 0.3, :V => 0.2)
            om = OmegaMatrix(omegas)

            @test length(om.param_names) == 2
            @test :CL in om.param_names
            @test :V in om.param_names

            # Check variances on diagonal (exact indices depend on sorting)
            cl_idx = findfirst(==(:CL), om.param_names)
            v_idx = findfirst(==(:V), om.param_names)
            @test om.matrix[cl_idx, cl_idx] ≈ 0.09 atol=1e-10  # 0.3^2
            @test om.matrix[v_idx, v_idx] ≈ 0.04 atol=1e-10    # 0.2^2

            # Off-diagonal should be zero
            @test om.matrix[1, 2] ≈ 0.0 atol=1e-10
            @test om.matrix[2, 1] ≈ 0.0 atol=1e-10

            # Cholesky should be lower triangular
            @test isapprox(om.cholesky_L, LowerTriangular(om.cholesky_L), atol=1e-10)
        end

        @testset "Full OmegaMatrix with Correlations" begin
            param_names = [:CL, :V]
            var_CL = 0.09  # SD = 0.3
            var_V = 0.04   # SD = 0.2
            corr = 0.5
            cov = corr * sqrt(var_CL) * sqrt(var_V)

            matrix = [var_CL cov; cov var_V]
            om = OmegaMatrix(param_names, matrix)

            @test om.param_names == [:CL, :V]
            @test isapprox(om.matrix[1, 1], var_CL, atol=1e-10)
            @test isapprox(om.matrix[2, 2], var_V, atol=1e-10)
            @test isapprox(om.matrix[1, 2], cov, atol=1e-10)
            @test isapprox(om.matrix[2, 1], cov, atol=1e-10)

            # Should be positive definite
            @test is_positive_definite(om.matrix)

            # Cholesky should satisfy L*L' = matrix
            reconstructed = om.cholesky_L * om.cholesky_L'
            @test isapprox(reconstructed, om.matrix, atol=1e-10)
        end

        @testset "Enforced Positive Definiteness" begin
            param_names = [:CL, :V]
            # Nearly singular matrix
            matrix = [1.0 0.999; 0.999 1.0]
            om = OmegaMatrix(param_names, matrix)

            # Should be positive definite after construction
            @test is_positive_definite(om.matrix)

            # Cholesky should exist
            @test om.cholesky_L !== nothing
        end

        @testset "get_diagonal_omegas" begin
            param_names = [:CL, :V]
            var_CL = 0.09
            var_V = 0.04
            matrix = [var_CL 0.0; 0.0 var_V]
            om = OmegaMatrix(param_names, matrix)

            diag_omegas = get_diagonal_omegas(om)

            @test diag_omegas[:CL] ≈ 0.3 atol=1e-10
            @test diag_omegas[:V] ≈ 0.2 atol=1e-10
        end

        @testset "get_correlation_matrix" begin
            param_names = [:CL, :V]
            var_CL = 0.09
            var_V = 0.04
            target_corr = 0.5
            cov = target_corr * sqrt(var_CL) * sqrt(var_V)
            matrix = [var_CL cov; cov var_V]

            om = OmegaMatrix(param_names, matrix)
            corr_mat = get_correlation_matrix(om)

            @test isapprox(corr_mat[1, 1], 1.0, atol=1e-10)
            @test isapprox(corr_mat[2, 2], 1.0, atol=1e-10)
            @test isapprox(corr_mat[1, 2], target_corr, atol=1e-10)
            @test isapprox(corr_mat[2, 1], target_corr, atol=1e-10)
        end
    end

    @testset "IIVSpec Construction" begin
        @testset "IIVSpec with Diagonal (Backward Compatible)" begin
            omegas = Dict(:CL => 0.3, :V => 0.2)
            iiv = IIVSpec(LogNormalIIV(), omegas, UInt64(12345), 10)

            @test iiv.kind isa LogNormalIIV
            @test iiv.omegas == omegas
            @test iiv.omega_matrix === nothing
            @test !has_correlations(iiv)
            @test iiv.seed == UInt64(12345)
            @test iiv.n == 10
        end

        @testset "IIVSpec with Full Covariance" begin
            param_names = [:CL, :V]
            matrix = [0.09 0.03; 0.03 0.04]
            om = OmegaMatrix(param_names, matrix)

            iiv = IIVSpec(LogNormalIIV(), om, UInt64(12345), 10)

            @test iiv.kind isa LogNormalIIV
            @test has_correlations(iiv)
            @test iiv.omega_matrix !== nothing
            @test iiv.omega_matrix.param_names == param_names

            # Diagonal omegas should be derived
            @test isapprox(iiv.omegas[:CL], 0.3, atol=1e-10)
            @test isapprox(iiv.omegas[:V], 0.2, atol=1e-10)
        end

        @testset "get_omega_matrix" begin
            # From diagonal spec
            omegas_diag = Dict(:CL => 0.3, :V => 0.2)
            iiv_diag = IIVSpec(LogNormalIIV(), omegas_diag, UInt64(12345), 10)
            om_diag = get_omega_matrix(iiv_diag)
            @test length(om_diag.param_names) == 2

            # From full covariance spec
            param_names = [:CL, :V]
            matrix = [0.09 0.03; 0.03 0.04]
            om_full = OmegaMatrix(param_names, matrix)
            iiv_full = IIVSpec(LogNormalIIV(), om_full, UInt64(12345), 10)
            om_returned = get_omega_matrix(iiv_full)
            @test om_returned === om_full
        end
    end

    @testset "Correlated Eta Sampling" begin
        @testset "Sampling Produces Correct Covariance Structure" begin
            # Create OmegaMatrix with known correlation
            param_names = [:CL, :V]
            var_CL = 0.09
            var_V = 0.04
            target_corr = 0.6
            cov = target_corr * sqrt(var_CL) * sqrt(var_V)
            matrix = [var_CL cov; cov var_V]

            om = OmegaMatrix(param_names, matrix)
            iiv = IIVSpec(LogNormalIIV(), om, UInt64(12345), 1000)

            # Sample many etas
            n_samples = 10000
            etas_CL = zeros(n_samples)
            etas_V = zeros(n_samples)

            rng = StableRNG(12345)
            for i in 1:n_samples
                eta_vec, _ = NeoPKPD.sample_eta_vector(om, rng)
                etas_CL[i] = eta_vec[1]
                etas_V[i] = eta_vec[2]
            end

            # Check sample means are approximately 0
            @test abs(mean(etas_CL)) < 0.05
            @test abs(mean(etas_V)) < 0.05

            # Check sample variances match target
            @test isapprox(var(etas_CL), var_CL, rtol=0.1)
            @test isapprox(var(etas_V), var_V, rtol=0.1)

            # Check sample correlation matches target
            sample_corr = cor(etas_CL, etas_V)
            @test isapprox(sample_corr, target_corr, atol=0.1)
        end

        @testset "Independent Sampling with Diagonal" begin
            omegas = Dict(:CL => 0.3, :V => 0.2)
            om = OmegaMatrix(omegas)

            n_samples = 10000
            etas_CL = zeros(n_samples)
            etas_V = zeros(n_samples)

            rng = StableRNG(12345)
            for i in 1:n_samples
                eta_vec, param_names = NeoPKPD.sample_eta_vector(om, rng)
                cl_idx = findfirst(==(:CL), param_names)
                v_idx = findfirst(==(:V), param_names)
                etas_CL[i] = eta_vec[cl_idx === nothing ? 1 : cl_idx]
                etas_V[i] = eta_vec[v_idx === nothing ? 2 : v_idx]
            end

            # Correlation should be approximately 0 for diagonal
            sample_corr = cor(etas_CL, etas_V)
            @test abs(sample_corr) < 0.05
        end

        @testset "Reproducibility with Seed" begin
            param_names = [:CL, :V]
            matrix = [0.09 0.03; 0.03 0.04]
            om = OmegaMatrix(param_names, matrix)

            # Sample with same seed twice
            rng1 = StableRNG(42)
            eta1, _ = NeoPKPD.sample_eta_vector(om, rng1)

            rng2 = StableRNG(42)
            eta2, _ = NeoPKPD.sample_eta_vector(om, rng2)

            @test eta1 ≈ eta2
        end
    end

    @testset "Positive Definiteness Enforcement" begin
        @testset "ensure_positive_definite_omega" begin
            # Already positive definite
            pd_matrix = [1.0 0.5; 0.5 1.0]
            result_pd = ensure_positive_definite_omega(pd_matrix)
            @test is_positive_definite(result_pd)
            @test isapprox(result_pd, pd_matrix, atol=1e-10)

            # Not positive definite
            not_pd = [1.0 1.5; 1.5 1.0]
            result_not_pd = ensure_positive_definite_omega(not_pd)
            @test is_positive_definite(result_not_pd)
        end

        @testset "OmegaMatrix Handles Non-PD Input" begin
            param_names = [:CL, :V]
            # Matrix with eigenvalue near zero
            matrix = [1.0 0.999; 0.999 1.0]

            # Should not throw
            om = OmegaMatrix(param_names, matrix)
            @test is_positive_definite(om.matrix)
        end
    end

    @testset "Serialization Round-Trip" begin
        @testset "Diagonal IIVSpec" begin
            omegas = Dict(:CL => 0.3, :V => 0.2)
            iiv_orig = IIVSpec(LogNormalIIV(), omegas, UInt64(12345), 10)

            # Serialize
            serialized = NeoPKPD._serialize_iiv(iiv_orig)

            @test haskey(serialized, "kind")
            @test haskey(serialized, "omegas")
            @test haskey(serialized, "seed")
            @test haskey(serialized, "n")
            @test !haskey(serialized, "omega_matrix")

            # Deserialize
            iiv_parsed = NeoPKPD._parse_iiv(serialized)

            @test iiv_parsed.omegas[:CL] ≈ iiv_orig.omegas[:CL]
            @test iiv_parsed.omegas[:V] ≈ iiv_orig.omegas[:V]
            @test !has_correlations(iiv_parsed)
        end

        @testset "Full Covariance IIVSpec" begin
            param_names = [:CL, :V]
            matrix = [0.09 0.03; 0.03 0.04]
            om = OmegaMatrix(param_names, matrix)
            iiv_orig = IIVSpec(LogNormalIIV(), om, UInt64(12345), 10)

            # Serialize
            serialized = NeoPKPD._serialize_iiv(iiv_orig)

            @test haskey(serialized, "omega_matrix")
            @test serialized["omega_matrix"]["param_names"] == ["CL", "V"]
            @test serialized["omega_matrix"]["n_params"] == 2

            # Deserialize
            iiv_parsed = NeoPKPD._parse_iiv(serialized)

            @test has_correlations(iiv_parsed)
            @test iiv_parsed.omega_matrix.param_names == param_names
            @test isapprox(iiv_parsed.omega_matrix.matrix, matrix, atol=1e-10)
        end
    end

    @testset "Three Parameter Covariance" begin
        # Test with 3 correlated random effects
        param_names = [:CL, :V, :Ka]

        # Create covariance matrix with various correlations
        # Variances: CL=0.09 (30% CV), V=0.04 (20% CV), Ka=0.25 (50% CV)
        var_CL = 0.09
        var_V = 0.04
        var_Ka = 0.25

        # Correlations: CL-V=0.5, CL-Ka=0.3, V-Ka=-0.2
        corr_CL_V = 0.5
        corr_CL_Ka = 0.3
        corr_V_Ka = -0.2

        cov_CL_V = corr_CL_V * sqrt(var_CL) * sqrt(var_V)
        cov_CL_Ka = corr_CL_Ka * sqrt(var_CL) * sqrt(var_Ka)
        cov_V_Ka = corr_V_Ka * sqrt(var_V) * sqrt(var_Ka)

        matrix = [
            var_CL    cov_CL_V  cov_CL_Ka;
            cov_CL_V  var_V     cov_V_Ka;
            cov_CL_Ka cov_V_Ka  var_Ka
        ]

        om = OmegaMatrix(param_names, matrix)

        @test length(om.param_names) == 3
        @test is_positive_definite(om.matrix)

        # Check correlation matrix
        corr_mat = get_correlation_matrix(om)
        @test isapprox(corr_mat[1, 2], corr_CL_V, atol=0.01)
        @test isapprox(corr_mat[1, 3], corr_CL_Ka, atol=0.01)
        @test isapprox(corr_mat[2, 3], corr_V_Ka, atol=0.01)

        # Sample and verify structure
        n_samples = 5000
        etas = zeros(n_samples, 3)
        rng = StableRNG(42)

        for i in 1:n_samples
            eta_vec, _ = NeoPKPD.sample_eta_vector(om, rng)
            etas[i, :] = eta_vec
        end

        # Check variances
        @test isapprox(var(etas[:, 1]), var_CL, rtol=0.15)
        @test isapprox(var(etas[:, 2]), var_V, rtol=0.15)
        @test isapprox(var(etas[:, 3]), var_Ka, rtol=0.15)

        # Check correlations
        sample_corr_CL_V = cor(etas[:, 1], etas[:, 2])
        sample_corr_CL_Ka = cor(etas[:, 1], etas[:, 3])
        sample_corr_V_Ka = cor(etas[:, 2], etas[:, 3])

        @test isapprox(sample_corr_CL_V, corr_CL_V, atol=0.1)
        @test isapprox(sample_corr_CL_Ka, corr_CL_Ka, atol=0.1)
        @test isapprox(sample_corr_V_Ka, corr_V_Ka, atol=0.1)
    end

end
