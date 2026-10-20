using Test
using NeoPKPD
using LinearAlgebra

@testset "Time-Varying Covariates" begin

    @testset "Interpolation Types" begin
        @test StepInterpolation() isa CovariateInterpolation
        @test LinearInterpolation() isa CovariateInterpolation
        @test LOCFInterpolation() isa CovariateInterpolation
        @test NOCBInterpolation() isa CovariateInterpolation
        @test SplineInterpolation() isa CovariateInterpolation
    end

    @testset "Transform Types" begin
        @test NoTransform() isa CovariateTransform
        @test LogTransform(70.0) isa CovariateTransform
        @test PowerTransform(70.0, 0.75) isa CovariateTransform
        @test NormalizedTransform(100.0, 50.0) isa CovariateTransform
    end

    @testset "CovariateTimeSeries" begin
        @testset "Basic construction" begin
            series = CovariateTimeSeries(
                :WT,
                [0.0, 24.0, 48.0],
                [70.0, 71.0, 72.0],
                StepInterpolation()
            )
            @test series.name == :WT
            @test length(series.times) == 3
            @test length(series.values) == 3
            @test series.interpolation isa StepInterpolation
        end

        @testset "With units" begin
            series = CovariateTimeSeries(
                :WT,
                [0.0, 24.0],
                [70.0, 72.0],
                LinearInterpolation();
                units="kg"
            )
            @test series.units == "kg"
        end

        @testset "Validation" begin
            # Must have matching lengths
            @test_throws AssertionError CovariateTimeSeries(
                :WT,
                [0.0, 24.0],
                [70.0],  # Mismatched length
                StepInterpolation()
            )

            # Must be sorted
            @test_throws AssertionError CovariateTimeSeries(
                :WT,
                [24.0, 0.0],  # Not sorted
                [70.0, 71.0],
                StepInterpolation()
            )
        end
    end

    @testset "Step Interpolation" begin
        series = CovariateTimeSeries(
            :WT,
            [0.0, 24.0, 48.0],
            [70.0, 72.0, 68.0],
            StepInterpolation()
        )

        # Before first measurement - use first value
        @test get_covariate_at_time(series, -1.0) == 70.0

        # At measurement times
        @test get_covariate_at_time(series, 0.0) == 70.0
        @test get_covariate_at_time(series, 24.0) == 72.0
        @test get_covariate_at_time(series, 48.0) == 68.0

        # Between measurements - use last value (LOCF)
        @test get_covariate_at_time(series, 12.0) == 70.0
        @test get_covariate_at_time(series, 36.0) == 72.0

        # After last measurement - use last value
        @test get_covariate_at_time(series, 100.0) == 68.0
    end

    @testset "Linear Interpolation" begin
        series = CovariateTimeSeries(
            :WT,
            [0.0, 24.0, 48.0],
            [70.0, 74.0, 72.0],
            LinearInterpolation()
        )

        # At measurement times
        @test get_covariate_at_time(series, 0.0) == 70.0
        @test get_covariate_at_time(series, 24.0) == 74.0
        @test get_covariate_at_time(series, 48.0) == 72.0

        # Midpoints - linear interpolation
        @test isapprox(get_covariate_at_time(series, 12.0), 72.0, atol=1e-10)  # (70 + 74) / 2
        @test isapprox(get_covariate_at_time(series, 36.0), 73.0, atol=1e-10)  # (74 + 72) / 2

        # Quarter points
        @test isapprox(get_covariate_at_time(series, 6.0), 71.0, atol=1e-10)  # 70 + 0.25 * (74 - 70)
    end

    @testset "NOCB Interpolation" begin
        series = CovariateTimeSeries(
            :WT,
            [0.0, 24.0, 48.0],
            [70.0, 72.0, 68.0],
            NOCBInterpolation()
        )

        # Between measurements - use next value
        @test get_covariate_at_time(series, 12.0) == 72.0
        @test get_covariate_at_time(series, 36.0) == 68.0
    end

    @testset "SubjectTimeVaryingCovariates" begin
        @testset "Construction" begin
            wt_series = CovariateTimeSeries(
                :WT,
                [0.0, 168.0],
                [70.0, 68.0],
                LinearInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series),
                baseline_covariates=Dict(:SEX => 1.0, :AGE => 45.0)
            )

            @test subject_cov.subject_id == "SUBJ001"
            @test haskey(subject_cov.covariates, :WT)
            @test haskey(subject_cov.baseline_covariates, :SEX)
            @test subject_cov.baseline_covariates[:AGE] == 45.0
        end

        @testset "Get covariates at times" begin
            wt_series = CovariateTimeSeries(
                :WT,
                [0.0, 24.0],
                [70.0, 72.0],
                LinearInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series),
                baseline_covariates=Dict(:SEX => 1.0)
            )

            # Time-varying covariate
            wt_values = get_covariates_at_times(subject_cov, :WT, [0.0, 12.0, 24.0])
            @test wt_values[1] == 70.0
            @test isapprox(wt_values[2], 71.0, atol=1e-10)
            @test wt_values[3] == 72.0

            # Baseline covariate - same at all times
            sex_values = get_covariates_at_times(subject_cov, :SEX, [0.0, 12.0, 24.0])
            @test all(sex_values .== 1.0)
        end
    end

    @testset "Transform Application" begin
        @testset "NoTransform" begin
            @test NeoPKPD.apply_transform(70.0, NoTransform()) == 70.0
        end

        @testset "LogTransform" begin
            transform = LogTransform(70.0)
            @test isapprox(NeoPKPD.apply_transform(70.0, transform), 0.0, atol=1e-10)
            @test isapprox(NeoPKPD.apply_transform(140.0, transform), log(2.0), atol=1e-10)
        end

        @testset "PowerTransform" begin
            transform = PowerTransform(70.0, 0.75)
            @test isapprox(NeoPKPD.apply_transform(70.0, transform), 1.0, atol=1e-10)
            @test isapprox(NeoPKPD.apply_transform(140.0, transform), 2.0^0.75, atol=1e-10)
        end

        @testset "NormalizedTransform" begin
            transform = NormalizedTransform(100.0, 50.0)
            @test isapprox(NeoPKPD.apply_transform(100.0, transform), 0.0, atol=1e-10)
            @test isapprox(NeoPKPD.apply_transform(150.0, transform), 1.0, atol=1e-10)
            @test isapprox(NeoPKPD.apply_transform(50.0, transform), -1.0, atol=1e-10)
        end
    end

    @testset "TimeVaryingCovariateEffect" begin
        @testset "Basic construction" begin
            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:WT,
                transform=PowerTransform(70.0, 0.75)
            )

            @test effect.parameter == :CL
            @test effect.parameter_index == 1
            @test effect.covariate == :WT
            @test effect.transform isa PowerTransform
            @test effect.fixed_effect == 1.0  # Default
        end

        @testset "With theta_index" begin
            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:CRCL,
                transform=NormalizedTransform(100.0, 50.0),
                theta_index=3
            )

            @test effect.theta_index == 3
            @test effect.fixed_effect === nothing
        end
    end

    @testset "TimeVaryingCovariateConfig" begin
        effect1 = TimeVaryingCovariateEffect(
            parameter=:CL,
            parameter_index=1,
            covariate=:WT,
            transform=PowerTransform(70.0, 0.75)
        )

        effect2 = TimeVaryingCovariateEffect(
            parameter=:V,
            parameter_index=2,
            covariate=:WT,
            transform=PowerTransform(70.0, 1.0)
        )

        config = TimeVaryingCovariateConfig(
            effects=[effect1, effect2],
            center_at_baseline=true,
            handle_missing=:locf
        )

        @test length(config.effects) == 2
        @test config.center_at_baseline == true
        @test config.handle_missing == :locf
    end

    @testset "Apply Time-Varying Covariates" begin
        @testset "Power transform (allometric)" begin
            # Weight effect on CL: CL * (WT/70)^0.75
            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:WT,
                transform=PowerTransform(70.0, 0.75)
            )

            wt_series = CovariateTimeSeries(
                :WT,
                [0.0, 24.0],
                [70.0, 140.0],  # Weight doubles
                StepInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series)
            )

            theta = [10.0, 50.0]  # CL=10, V=50

            # At t=0, WT=70: CL should be unchanged
            theta_t0 = apply_time_varying_covariates(theta, 0.0, [effect], subject_cov)
            @test isapprox(theta_t0[1], 10.0, atol=1e-10)  # CL * (70/70)^0.75 = CL

            # At t=24, WT=140: CL should increase by 2^0.75
            theta_t24 = apply_time_varying_covariates(theta, 24.0, [effect], subject_cov)
            @test isapprox(theta_t24[1], 10.0 * 2.0^0.75, atol=1e-10)

            # V should be unchanged
            @test theta_t0[2] == 50.0
            @test theta_t24[2] == 50.0
        end

        @testset "Normalized transform" begin
            # CRCL effect: CL * (1 + θ_CRCL * (CRCL - 100)/50)
            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:CRCL,
                transform=NormalizedTransform(100.0, 50.0),
                theta_index=3
            )

            crcl_series = CovariateTimeSeries(
                :CRCL,
                [0.0, 24.0],
                [100.0, 150.0],
                LinearInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:CRCL => crcl_series)
            )

            theta = [10.0, 50.0, 0.5]  # CL=10, V=50, θ_CRCL=0.5

            # At t=0, CRCL=100: CL should be unchanged
            theta_t0 = apply_time_varying_covariates(theta, 0.0, [effect], subject_cov)
            @test isapprox(theta_t0[1], 10.0, atol=1e-10)  # CL * (1 + 0.5 * 0) = CL

            # At t=24, CRCL=150: CL * (1 + 0.5 * 1) = CL * 1.5
            theta_t24 = apply_time_varying_covariates(theta, 24.0, [effect], subject_cov)
            @test isapprox(theta_t24[1], 10.0 * 1.5, atol=1e-10)
        end

        @testset "Multiple effects" begin
            effects = [
                TimeVaryingCovariateEffect(
                    parameter=:CL,
                    parameter_index=1,
                    covariate=:WT,
                    transform=PowerTransform(70.0, 0.75)
                ),
                TimeVaryingCovariateEffect(
                    parameter=:V,
                    parameter_index=2,
                    covariate=:WT,
                    transform=PowerTransform(70.0, 1.0)
                )
            ]

            wt_series = CovariateTimeSeries(
                :WT,
                [0.0],
                [140.0],  # 2x reference weight
                StepInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series)
            )

            theta = [10.0, 50.0]

            theta_mod = apply_time_varying_covariates(theta, 0.0, effects, subject_cov)

            # CL scaled by 2^0.75
            @test isapprox(theta_mod[1], 10.0 * 2.0^0.75, atol=1e-10)

            # V scaled by 2^1.0
            @test isapprox(theta_mod[2], 50.0 * 2.0, atol=1e-10)
        end
    end

    @testset "Convenience Constructors" begin
        @testset "weight_effect_on_cl" begin
            effect = weight_effect_on_cl(reference=70.0, power=0.75)
            @test effect.parameter == :CL
            @test effect.covariate == :WT
            @test effect.transform isa PowerTransform
            @test effect.transform.reference == 70.0
            @test effect.transform.power == 0.75
        end

        @testset "weight_effect_on_v" begin
            effect = weight_effect_on_v(reference=70.0, power=1.0)
            @test effect.parameter == :V
            @test effect.covariate == :WT
            @test effect.transform.power == 1.0
        end

        @testset "crcl_effect_on_cl" begin
            effect = crcl_effect_on_cl(reference=100.0, scale=50.0, theta_index=3)
            @test effect.parameter == :CL
            @test effect.covariate == :CRCL
            @test effect.theta_index == 3
            @test effect.transform isa NormalizedTransform
        end

        @testset "albumin_effect_on_cl" begin
            effect = albumin_effect_on_cl(reference=4.0)
            @test effect.covariate == :ALB
            @test effect.transform isa PowerTransform
            @test effect.transform.reference == 4.0
        end
    end

    @testset "Predictions with Time-Varying Covariates" begin
        @testset "One-compartment IV bolus" begin
            # Model spec
            model_params = OneCompIVBolusParams(10.0, 50.0)
            doses = [DoseEvent(0.0, 100.0)]
            model_spec = ModelSpec(OneCompIVBolus(), "pk_iv", model_params, doses)

            # Weight effect on CL
            effects = [
                TimeVaryingCovariateEffect(
                    parameter=:CL,
                    parameter_index=1,
                    covariate=:WT,
                    transform=PowerTransform(70.0, 0.75)
                )
            ]

            # Time-varying weight (decreasing over time - e.g., cancer patient)
            wt_series = CovariateTimeSeries(
                :WT,
                [0.0, 48.0, 96.0],
                [70.0, 65.0, 60.0],
                LinearInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series)
            )

            theta = [10.0, 50.0]
            eta = [0.0, 0.0]
            times = [0.5, 12.0, 48.0, 96.0]

            pred = compute_predictions_with_tv_covariates(
                theta, eta, times, doses, model_spec, effects, subject_cov
            )

            @test length(pred) == 4
            @test all(isfinite.(pred))
            @test all(pred .> 0)

            # Later predictions should be higher (slower elimination due to lower weight)
            # Because CL decreases as WT decreases
        end
    end

    @testset "Validation" begin
        @testset "Covariate coverage check" begin
            wt_series = CovariateTimeSeries(
                :WT,
                [24.0, 48.0],  # Starts at t=24
                [70.0, 72.0],
                StepInterpolation()
            )

            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001";
                covariates=Dict(:WT => wt_series)
            )

            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:WT,
                transform=PowerTransform(70.0, 0.75)
            )

            # Observations starting before covariate data
            obs_times = [0.0, 12.0, 24.0, 48.0, 72.0]

            warnings = validate_covariate_coverage(subject_cov, obs_times, [effect])

            @test length(warnings) >= 1
            @test any(contains(w, "not available before") for w in warnings)
        end

        @testset "Missing covariate detection" begin
            subject_cov = SubjectTimeVaryingCovariates(
                "SUBJ001"  # No covariates
            )

            effect = TimeVaryingCovariateEffect(
                parameter=:CL,
                parameter_index=1,
                covariate=:CRCL,
                transform=NormalizedTransform(100.0, 50.0)
            )

            warnings = validate_covariate_coverage(subject_cov, [0.0, 24.0], [effect])

            @test length(warnings) >= 1
            @test any(contains(w, "not found") for w in warnings)
        end
    end

    @testset "Extract from SubjectData" begin
        @testset "Scalar covariates" begin
            subject = SubjectData(
                "SUBJ001",
                [0.0, 1.0, 2.0],
                [10.0, 8.0, 6.0],
                [DoseEvent(0.0, 100.0)];
                covariates=Dict{Symbol,Any}(:SEX => 1, :AGE => 45)
            )

            subject_cov = extract_time_varying_covariates(subject)

            @test subject_cov.subject_id == "SUBJ001"
            @test haskey(subject_cov.baseline_covariates, :SEX)
            @test haskey(subject_cov.baseline_covariates, :AGE)
            @test subject_cov.baseline_covariates[:SEX] == 1.0
            @test subject_cov.baseline_covariates[:AGE] == 45.0
        end

        @testset "Vector covariates (time-varying at obs times)" begin
            subject = SubjectData(
                "SUBJ001",
                [0.0, 24.0, 48.0],
                [10.0, 8.0, 6.0],
                [DoseEvent(0.0, 100.0)];
                covariates=Dict{Symbol,Any}(:WT => [70.0, 71.0, 72.0])
            )

            subject_cov = extract_time_varying_covariates(subject)

            @test haskey(subject_cov.covariates, :WT)
            wt_series = subject_cov.covariates[:WT]
            @test wt_series.times == [0.0, 24.0, 48.0]
            @test wt_series.values == [70.0, 71.0, 72.0]
        end
    end

    @testset "Pretty Printing" begin
        @testset "TimeVaryingCovariateEffect show" begin
            effect = weight_effect_on_cl()
            io = IOBuffer()
            show(io, effect)
            output = String(take!(io))
            @test contains(output, "TimeVaryingCovariateEffect")
            @test contains(output, "CL")
            @test contains(output, "WT")
        end

        @testset "CovariateTimeSeries show" begin
            series = CovariateTimeSeries(
                :WT,
                [0.0, 24.0, 48.0],
                [70.0, 72.0, 74.0],
                LinearInterpolation()
            )
            io = IOBuffer()
            show(io, series)
            output = String(take!(io))
            @test contains(output, "CovariateTimeSeries")
            @test contains(output, "WT")
            @test contains(output, "3 measurements")
        end
    end

end
