@testset "Trial Module" begin

    @testset "Study Designs" begin
        @testset "Parallel Design" begin
            # Basic parallel design
            design = parallel_design(2)
            @test design isa ParallelDesign
            @test design.n_arms == 2
            @test length(design.randomization_ratio) == 2
            @test sum(design.randomization_ratio) ≈ 1.0

            # 3-arm with custom ratio
            design3 = parallel_design(3, randomization_ratio=[0.5, 0.25, 0.25])
            @test design3.n_arms == 3
            @test design3.randomization_ratio ≈ [0.5, 0.25, 0.25]
        end

        @testset "Crossover Design" begin
            # 2x2 crossover
            design = crossover_2x2()
            @test design isa CrossoverDesign
            @test design.n_periods == 2
            @test design.n_sequences == 2
            @test design.washout_duration == 7.0

            # Williams design
            williams = williams_design(4)
            @test williams.n_periods == 4
            @test length(williams.sequence_assignments) == 8
        end

        @testset "Dose Escalation" begin
            # 3+3 design
            dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
            design = dose_escalation_3plus3(dose_levels)
            @test design isa DoseEscalationDesign
            @test design.starting_dose == 10.0
            @test design.cohort_size == 3
            @test design.escalation_rule isa ThreePlusThree

            # mTPI design
            mtpi = dose_escalation_mtpi(dose_levels, target_dlt_rate=0.30)
            @test mtpi.escalation_rule isa mTPI
            @test mtpi.escalation_rule.target_dlt_rate == 0.30

            # CRM design
            crm = dose_escalation_crm(dose_levels)
            @test crm.escalation_rule isa CRM
        end

        @testset "Bioequivalence Design" begin
            design = bioequivalence_design()
            @test design isa BioequivalenceDesign
            @test design.bioequivalence_limits == (0.80, 1.25)
            @test :cmax in design.parameters
            @test :auc_0_inf in design.parameters
        end

        @testset "Adaptive Design" begin
            base = parallel_design(2)
            adaptive = adaptive_design(base, interim_analyses=[0.5])
            @test adaptive isa AdaptiveDesign
            @test adaptive.interim_analyses == [0.5]
            @test adaptive.alpha_spending == :obrien_fleming
        end

        @testset "Design Descriptions" begin
            @test contains(get_design_description(parallel_design(2)), "parallel")
            @test contains(get_design_description(crossover_2x2()), "crossover")
            @test contains(get_design_description(bioequivalence_design()), "Bioequivalence")
        end
    end

    @testset "Dosing Regimens" begin
        @testset "Standard Regimens" begin
            qd = dosing_qd(100.0, 7)
            @test qd.frequency isa QD
            @test qd.dose_amount == 100.0
            @test qd.duration_days == 7

            bid = dosing_bid(50.0, 14)
            @test bid.frequency isa BID
            @test length(bid.dose_times) == 2

            tid = dosing_tid(25.0, 7)
            @test tid.frequency isa TID
            @test length(tid.dose_times) == 3

            # With loading dose
            with_loading = dosing_qd(100.0, 7, loading_dose=200.0)
            @test with_loading.loading_dose == 200.0
        end

        @testset "Dose Event Times" begin
            regimen = dosing_bid(50.0, 3)
            times = dose_event_times(regimen)
            @test length(times) == 6  # 3 days * 2 doses/day
            @test times[1] == 8.0  # First dose at 8 AM
        end

        @testset "Titration Regimen" begin
            tit = titration_regimen(25.0, 100.0, 4, 7)
            @test tit isa TitrationRegimen
            @test length(tit.steps) == 4
            @test tit.steps[1].dose == 25.0
            @test tit.steps[end].dose == 100.0

            duration = total_regimen_duration(tit)
            @test duration == 28  # 4 steps * 7 days
        end

        @testset "Generate Doses" begin
            regimen = dosing_qd(100.0, 5)
            doses = generate_doses(regimen)
            @test length(doses) == 5
            @test all(d == 100.0 for d in doses)

            # With loading dose
            regimen_load = dosing_qd(100.0, 5, loading_dose=200.0)
            doses_load = generate_doses(regimen_load)
            @test doses_load[1] == 200.0
            @test doses_load[2] == 100.0
        end
    end

    @testset "Virtual Population" begin
        @testset "Demographic Specs" begin
            demo = default_demographic_spec()
            @test demo.age_mean == 35.0
            @test demo.weight_mean == 75.0

            hv = healthy_volunteer_spec()
            @test hv.age_range == (18.0, 45.0)
        end

        @testset "Patient Populations" begin
            demo, disease = patient_population_spec(:diabetes)
            @test demo.age_mean > 40.0  # Older population
            @test disease.name == :diabetes

            demo2, disease2 = patient_population_spec(:renal)
            @test disease2.name == :renal_impairment
        end

        @testset "Population Generation" begin
            spec = VirtualPopulationSpec(
                demographics = DemographicSpec(age_mean=40.0, age_sd=10.0),
                seed = UInt64(42)
            )
            pop = generate_virtual_population(spec, 50)

            @test length(pop) == 50
            @test all(p.id >= 1 for p in pop)
            @test all(18.0 <= p.age <= 75.0 for p in pop)
            @test all(p.sex in [:male, :female] for p in pop)
        end

        @testset "Population Summary" begin
            spec = VirtualPopulationSpec(seed = UInt64(42))
            pop = generate_virtual_population(spec, 100)
            summary = summarize_population(pop)

            @test summary[:n] == 100
            @test haskey(summary, :age_mean)
            @test haskey(summary, :weight_mean)
            @test haskey(summary, :female_proportion)
            @test 0.0 <= summary[:female_proportion] <= 1.0
        end
    end

    @testset "Trial Events" begin
        @testset "Dropout Simulation" begin
            spec = DropoutSpec(random_rate_per_day=0.01)
            dropouts = simulate_dropout(spec, 30.0, 100)

            @test length(dropouts) <= 100
            @test all(d.time_days <= 30.0 for d in dropouts)
            @test all(d.reason in [:random, :ae, :non_compliance] for d in dropouts)
        end

        @testset "Compliance" begin
            spec = ComplianceSpec(mean_compliance=0.85)
            compliance = apply_compliance(spec, 50, 28.0)

            @test length(compliance) == 50
            @test all(0.0 <= c <= 1.0 for c in compliance)
            @test 0.7 < sum(compliance)/50 < 1.0  # Average around 0.85
        end

        @testset "Survival Time" begin
            dropouts = [DropoutEvent(1, 10.0, :random), DropoutEvent(3, 15.0, :ae)]

            time1, comp1 = calculate_survival_time(dropouts, 1, 30.0)
            @test time1 == 10.0
            @test !comp1

            time2, comp2 = calculate_survival_time(dropouts, 2, 30.0)
            @test time2 == 30.0
            @test comp2
        end
    end

    @testset "Endpoints" begin
        @testset "PK Endpoint Extraction" begin
            result = Dict{String, Any}(
                "t" => [0.0, 1.0, 2.0, 4.0, 8.0],
                "observations" => Dict{String, Vector{Float64}}(
                    "conc" => [0.0, 10.0, 8.0, 5.0, 2.0]
                )
            )

            endpoint_cmax = PKEndpoint(:cmax, metric=:cmax)
            cmax = extract_pk_endpoint(result, endpoint_cmax)
            @test cmax == 10.0

            endpoint_tmax = PKEndpoint(:tmax, metric=:tmax)
            tmax = extract_pk_endpoint(result, endpoint_tmax)
            @test tmax == 1.0

            endpoint_auc = PKEndpoint(:auc, metric=:auc_0_t)
            auc = extract_pk_endpoint(result, endpoint_auc)
            @test auc > 0
        end

        @testset "Endpoint Analysis" begin
            values = [10.0, 12.0, 9.0, 11.0, 13.0, 8.0, 10.5, 11.5, 9.5, 10.0]
            endpoint = PKEndpoint(:test, metric=:cmax)
            stats = analyze_endpoint(values, endpoint)

            @test stats[:n] == 10
            @test 9.0 < stats[:mean] < 12.0
            @test stats[:min] == 8.0
            @test stats[:max] == 13.0
        end

        @testset "Arm Comparison" begin
            arm1 = [10.0, 11.0, 12.0, 9.0, 10.5, 11.5, 10.0, 9.5]
            arm2 = [14.0, 15.0, 16.0, 13.0, 14.5, 15.5, 14.0, 13.5]

            comparison = compare_arms(arm1, arm2, test=:ttest)
            @test comparison[:difference] > 0  # arm2 > arm1
            @test haskey(comparison, :t_statistic)
            @test haskey(comparison, :ci_lower)
            @test haskey(comparison, :ci_upper)
        end

        @testset "Responder Analysis" begin
            values = [0.3, 0.6, 0.8, 0.4, 0.9, 0.2, 0.7, 0.5, 0.85, 0.35]
            result = responder_analysis(values, 0.5, direction=:greater)

            @test result[:n] == 10
            @test result[:n_responders] == 6  # values >= 0.5
            @test result[:response_rate] == 0.6
        end
    end

    @testset "Power Analysis" begin
        @testset "Analytical Power" begin
            # Standard two-sample t-test scenario
            power = estimate_power_analytical(50, 0.5, 1.0)
            @test 0.0 < power < 1.0

            # Larger effect = higher power
            power_large = estimate_power_analytical(50, 1.0, 1.0)
            @test power_large > power

            # Larger N = higher power
            power_large_n = estimate_power_analytical(100, 0.5, 1.0)
            @test power_large_n > power
        end

        @testset "Sample Size Estimation" begin
            result = estimate_sample_size(0.80, 0.5, 1.0)
            @test result isa SampleSizeResult
            @test result.n_per_arm > 0
            @test result.achieved_power >= result.target_power
        end

        @testset "Alpha Spending" begin
            # O'Brien-Fleming - conservative early
            alpha_50 = alpha_spending_function(0.5, 0.05, :obrien_fleming)
            @test alpha_50 < 0.025  # Less than half of total alpha

            # Linear spending
            alpha_lin = alpha_spending_function(0.5, 0.05, :linear)
            @test alpha_lin ≈ 0.025 atol=0.001

            # At final analysis
            alpha_final = alpha_spending_function(1.0, 0.05, :obrien_fleming)
            @test alpha_final ≈ 0.05 atol=0.01
        end

        @testset "Incremental Alpha" begin
            fractions = [0.5, 1.0]
            incremental = incremental_alpha(fractions, 0.05, :obrien_fleming)

            @test length(incremental) == 2
            @test sum(incremental) ≈ 0.05 atol=0.01
        end
    end

    @testset "Paired vs Unpaired T-Test Validation" begin
        # Generate correlated paired data where paired test should be more powerful
        rng = StableRNG(12345)

        # Create baseline values with substantial variation
        baseline = randn(rng, 30) .* 10.0 .+ 100.0

        # Treatment effect of 5 units with small within-subject noise
        effect = 5.0
        treatment = baseline .+ effect .+ randn(rng, 30) .* 3.0

        # Paired t-test
        paired_result = compare_arms(baseline, treatment; test=:ttest, paired=true)

        # Unpaired t-test (same data, but ignoring pairing)
        unpaired_result = compare_arms(baseline, treatment; test=:ttest, paired=false)

        @test paired_result[:test_type] == :paired_ttest
        @test unpaired_result[:test_type] == :welch_ttest

        # Paired test should have a smaller (more significant) p-value
        # because it removes between-subject variability
        @test paired_result[:p_value] < unpaired_result[:p_value]

        # Effect estimate should be close to true effect
        @test abs(paired_result[:difference] - effect) < 2.0

        # Degrees of freedom: paired = n-1, Welch = different
        @test paired_result[:df] == 29  # n-1
    end

    @testset "Crossover Analysis Module" begin
        @testset "Bioequivalence Assessment" begin
            # Generate data simulating two formulations with GMR = 1.0
            rng = StableRNG(54321)
            n = 24

            # Reference (log-normal, mean ~100)
            log_ref = randn(rng, n) .* 0.2 .+ log(100.0)
            reference = exp.(log_ref)

            # Test formulation: similar to reference (GMR ~1.0)
            log_test = log_ref .+ randn(rng, n) .* 0.15  # Add within-subject variability
            test_values = exp.(log_test)

            be_result = assess_bioequivalence(reference, test_values)

            @test haskey(be_result, :geometric_mean_ratio)
            @test haskey(be_result, :ci_90_lower)
            @test haskey(be_result, :ci_90_upper)
            @test haskey(be_result, :is_bioequivalent)

            # With similar formulations, GMR should be close to 1.0
            @test 0.9 < be_result[:geometric_mean_ratio] < 1.1

            # Test formulation that is NOT bioequivalent
            test_high = reference .* 1.4  # 40% higher
            be_fail = assess_bioequivalence(reference, test_high)
            @test be_fail[:is_bioequivalent] == false
        end

        @testset "Within-Subject CV" begin
            rng = StableRNG(11111)
            n = 20

            # Create paired data with known within-subject CV
            v1 = exp.(randn(rng, n) .* 0.3 .+ log(100.0))
            v2 = exp.(log.(v1) .+ randn(rng, n) .* 0.2)  # ~20% within-subject CV

            cv = compute_within_subject_cv(v1, v2)

            @test !isnan(cv)
            @test cv > 0
            @test cv < 50.0  # Should be reasonable CV
        end

        @testset "Period Effect Test" begin
            design = crossover_2x2()

            # Data with no period effect
            period1 = [100.0, 105.0, 98.0, 102.0, 110.0, 95.0, 108.0, 97.0]
            period2 = [101.0, 104.0, 99.0, 103.0, 109.0, 96.0, 107.0, 98.0]

            result = test_period_effect(period1, period2, design)

            @test haskey(result, :mean_period1)
            @test haskey(result, :mean_period2)
            @test haskey(result, :n_subjects)
            @test result[:n_subjects] == 8

            # With similar means, should not detect significant period effect
            if haskey(result, :p_value)
                @test result[:p_value] > 0.1
            end
        end
    end

    @testset "Adaptive Trial Module" begin
        @testset "Alpha Spending Functions" begin
            alpha = 0.025

            # O'Brien-Fleming: very conservative early
            obf_30 = get_alpha_spending(:obrien_fleming, 0.3, alpha)
            obf_50 = get_alpha_spending(:obrien_fleming, 0.5, alpha)
            obf_100 = get_alpha_spending(:obrien_fleming, 1.0, alpha)

            @test obf_30 < obf_50 < obf_100
            @test obf_30 < 0.005  # Very little alpha spent early
            @test obf_100 ≈ alpha atol=0.005

            # Pocock: more aggressive early spending
            pocock_50 = get_alpha_spending(:pocock, 0.5, alpha)
            @test pocock_50 > obf_50  # Pocock spends more early

            # Linear spending
            linear_50 = get_alpha_spending(:linear, 0.5, alpha)
            @test linear_50 ≈ 0.5 * alpha atol=0.002
        end

        @testset "Efficacy Boundary" begin
            alpha = 0.025

            # At final analysis (info_frac = 1.0), boundary should be z_alpha
            boundary_final = get_efficacy_boundary(:obrien_fleming, 1.0, alpha)
            z_alpha = 1.96  # Approximately
            @test abs(boundary_final - z_alpha) < 0.2

            # Early interim should have higher boundary (more stringent)
            boundary_early = get_efficacy_boundary(:obrien_fleming, 0.5, alpha)
            @test boundary_early > boundary_final
        end

        @testset "Stopping Rules" begin
            # Efficacy stop: z_stat exceeds boundary
            @test check_efficacy_stop(3.5, 2.8) == true
            @test check_efficacy_stop(2.5, 2.8) == false

            # Futility stop: conditional power too low
            @test check_futility_stop(0.05, 0.10) == true
            @test check_futility_stop(0.20, 0.10) == false
        end

        @testset "Conditional Power" begin
            effect = 0.4  # Observed effect
            se = 0.15     # Standard error
            info_frac = 0.5
            target = 0.5   # Assumed true effect

            cp = compute_conditional_power(effect, se, info_frac, target)

            @test !isnan(cp)
            @test 0.0 <= cp <= 1.0

            # Stronger observed effect should give higher conditional power
            cp_strong = compute_conditional_power(0.6, se, info_frac, target)
            @test cp_strong > cp
        end

        @testset "InterimResult Structure" begin
            interim = InterimResult(
                1,      # analysis_number
                0.5,    # information_fraction
                50,     # n_enrolled
                0.3,    # primary_effect
                0.12,   # primary_se
                2.5,    # primary_z_stat
                0.012,  # p_value
                2.8,    # efficacy_boundary
                0.1,    # futility_boundary
                0.75,   # conditional_power
                false,  # stop_for_efficacy
                false   # stop_for_futility
            )

            @test interim.analysis_number == 1
            @test interim.information_fraction == 0.5
            @test interim.n_enrolled == 50
            @test interim.conditional_power == 0.75
        end
    end

    @testset "Trial Specification Validation" begin
        @testset "Valid Parallel Trial" begin
            # Create a simple valid spec using the correct API
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            regimen = dosing_qd(100.0, 7)
            arm1 = TreatmentArm("Control", model_spec, regimen; n_subjects=25)
            arm2 = TreatmentArm("Treatment", model_spec, regimen; n_subjects=25)

            design = parallel_design(2)
            trial = TrialSpec("Valid Trial", design, [arm1, arm2];
                              endpoints = EndpointSpec[PKEndpoint(:auc)])

            result = validate_trial_spec(trial)

            @test result.valid == true
            @test isempty(result.errors)
        end

        @testset "Invalid - Missing Model Spec" begin
            regimen = dosing_qd(100.0, 7)
            arm1 = TreatmentArm("Control", nothing, regimen; n_subjects=25)
            arm2 = TreatmentArm("Treatment", nothing, regimen; n_subjects=25)

            design = parallel_design(2)
            trial = TrialSpec("Invalid Trial", design, [arm1, arm2])

            result = validate_trial_spec(trial)

            @test result.valid == false
            @test any(contains(e, "pk_model_spec") for e in result.errors)
        end

        @testset "Invalid - Unimplemented Adaptive Feature" begin
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            regimen = dosing_qd(100.0, 7)
            arm = TreatmentArm("Treatment", model_spec, regimen; n_subjects=50)

            base = parallel_design(2)

            # Test platform trial (unimplemented)
            adaptive_platform = AdaptiveDesign(base;
                interim_analyses = [0.5],
                adaptation_rules = Dict{Symbol, Any}(
                    :platform_trial => true
                )
            )

            trial_platform = TrialSpec("Platform Trial", adaptive_platform, [arm, arm])
            result_platform = validate_trial_spec(trial_platform)

            @test result_platform.valid == false
            @test any(contains(e, "Platform/umbrella trials") for e in result_platform.errors)

            # Test basket trial (unimplemented)
            adaptive_basket = AdaptiveDesign(base;
                interim_analyses = [0.5],
                adaptation_rules = Dict{Symbol, Any}(
                    :basket_trial => true
                )
            )

            trial_basket = TrialSpec("Basket Trial", adaptive_basket, [arm, arm])
            result_basket = validate_trial_spec(trial_basket)

            @test result_basket.valid == false
            @test any(contains(e, "Basket trials") for e in result_basket.errors)
        end

        @testset "Strict Mode Warnings" begin
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            regimen = dosing_qd(100.0, 7)
            arm = TreatmentArm("Treatment", model_spec, regimen; n_subjects=25)

            design = parallel_design(2)
            # No endpoints defined - should generate warning
            trial = TrialSpec("Warning Trial", design, [arm, arm])

            # Non-strict: valid but with warnings
            result_normal = validate_trial_spec(trial)
            @test result_normal.valid == true
            @test !isempty(result_normal.warnings)

            # Strict: warnings become errors
            result_strict = validate_trial_spec(trial; strict = true)
            @test result_strict.valid == false
        end

        @testset "BioequivalenceDesign Validation" begin
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            regimen = dosing_qd(100.0, 1)
            arm1 = TreatmentArm("Reference", model_spec, regimen; n_subjects=12)
            arm2 = TreatmentArm("Test", model_spec, regimen; n_subjects=12)

            design = bioequivalence_design()

            # Valid BE trial with PK endpoints
            trial_valid = TrialSpec("BE Trial", design, [arm1, arm2];
                                    endpoints = EndpointSpec[PKEndpoint(:auc), PKEndpoint(:cmax)])
            result = validate_trial_spec(trial_valid)
            @test result.valid == true

            # Invalid: no PK endpoints
            trial_no_pk = TrialSpec("BE Trial No PK", design, [arm1, arm2];
                                    endpoints = EndpointSpec[PDEndpoint(:response)])
            result_no_pk = validate_trial_spec(trial_no_pk)
            @test result_no_pk.valid == false
            @test any(contains(e, "PKEndpoint") for e in result_no_pk.errors)
        end

        @testset "CRM Skeleton Validation" begin
            # Valid CRM skeleton
            dose_levels = [10.0, 25.0, 50.0, 100.0, 200.0]
            valid_skeleton = [0.05, 0.10, 0.20, 0.35, 0.50]
            valid_crm = dose_escalation_crm(dose_levels, skeleton=valid_skeleton)

            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            regimen = dosing_qd(10.0, 1)
            arm = TreatmentArm("Dose", model_spec, regimen; n_subjects=30)

            trial_valid = TrialSpec("CRM Trial", valid_crm, [arm])
            result_valid = validate_trial_spec(trial_valid)
            @test result_valid.valid == true

            # Invalid: skeleton not monotonic
            bad_skeleton = [0.05, 0.30, 0.20, 0.35, 0.50]  # Not monotonic
            bad_crm = dose_escalation_crm(dose_levels, skeleton=bad_skeleton)
            trial_bad = TrialSpec("Bad CRM", bad_crm, [arm])
            result_bad = validate_trial_spec(trial_bad)
            @test result_bad.valid == false
            @test any(contains(e, "monotonically increasing") for e in result_bad.errors)
        end
    end

    @testset "Real PK Simulation Integration" begin
        @testset "Subject Exposure Simulation" begin
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            # Create a virtual subject using convenience constructor
            subject = VirtualSubject(
                1,      # id
                35.0,   # age
                70.0,   # weight
                :male,  # sex
                :caucasian  # race
            )

            dose_events = [DoseEvent(0.0, 100.0)]
            observation_times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]

            exposure = simulate_subject_exposure(
                model_spec,
                subject,
                dose_events,
                observation_times;
                include_iiv = false  # Deterministic for testing
            )

            @test exposure isa SubjectExposure
            @test length(exposure.times) == length(observation_times)
            @test length(exposure.concentrations) == length(observation_times)

            # Concentrations should follow expected PK profile
            @test exposure.concentrations[1] ≈ 0.0 atol=0.01  # Before absorption
            @test maximum(exposure.concentrations) > 0  # Should have some peak
        end

        @testset "IIV Application" begin
            pk_model = OneCompOralFirstOrder()
            pk_params = OneCompOralFirstOrderParams(1.2, 5.0, 50.0)
            model_spec = ModelSpec(pk_model, "pk", pk_params, DoseEvent[])

            dose_events = [DoseEvent(0.0, 100.0)]
            observation_times = [0.0, 1.0, 4.0, 12.0, 24.0]

            # Create a subject without IIV using convenience constructor
            subject = VirtualSubject(1, 35.0, 70.0, :male, :caucasian)

            # Test with include_iiv=false (no IIV applied)
            exposure_no_iiv = simulate_subject_exposure(model_spec, subject, dose_events, observation_times;
                                                         include_iiv = false)

            @test !isempty(exposure_no_iiv.concentrations)
            @test maximum(exposure_no_iiv.concentrations) > 0

            # Test with include_iiv=true but no eta values (should gracefully handle)
            exposure_with_flag = simulate_subject_exposure(model_spec, subject, dose_events, observation_times;
                                                            include_iiv = true)

            # Without eta values, should produce same results
            @test exposure_no_iiv.concentrations ≈ exposure_with_flag.concentrations
        end
    end

end
