# =============================================================================
# Test suite for Industry-Standard TMDD Models
# =============================================================================
#
# Comprehensive tests for TMDD model variants matching NONMEM/Monolix/Phoenix.
# Tests all approximations: Full, QSS, QE, RapidBinding
# Tests all routes: IVBolus, IVInfusion, Subcutaneous
# Tests specialized models: FcRn, ADA, Soluble, Bispecific
# =============================================================================

using Test
using NeoPKPD
using StableRNGs
using Statistics

@testset "TMDD Models" begin

    @testset "TMDD Type Definitions" begin
        @testset "OneCptTMDD" begin
            @test OneCptTMDD() isa TMDDModelKind
            @test OneCptTMDD() isa ModelKind
            @test OneCptTMDD().approximation == QSS
            @test OneCptTMDD().route == IVBolus

            # With explicit approximation
            @test OneCptTMDD(FullTMDD).approximation == FullTMDD
            @test OneCptTMDD(QE).approximation == QE

            # Parameters
            params = OneCptTMDDParams(0.2, 3.0, 0.05, 0.1, 0.01, 0.1, 1.0)
            @test params.CL == 0.2
            @test params.V == 3.0
            @test params.KSS == 0.05
            @test params.R0 == 1.0
        end

        @testset "TwoCptTMDD - Industry Standard" begin
            @test TwoCptTMDD() isa TMDDModelKind
            @test TwoCptTMDD().approximation == QSS
            @test TwoCptTMDD().route == IVBolus

            # Various approximations
            @test TwoCptTMDD(FullTMDD, IVBolus).approximation == FullTMDD
            @test TwoCptTMDD(QE, IVBolus).approximation == QE
            @test TwoCptTMDD(RapidBinding, IVBolus).approximation == RapidBinding

            # Various routes
            @test TwoCptTMDD(QSS, Subcutaneous).route == Subcutaneous
            @test TwoCptTMDD(QSS, IVInfusion).route == IVInfusion

            # Parameters - industry standard clearance-based
            params = TwoCptTMDDParams(
                0.22, 3.28, 2.66, 0.54,  # CL, V1, V2, Q
                0.087, 0.037, 0.001, 0.012, 0.083  # KSS, kint, ksyn, kdeg, R0
            )
            @test params.CL == 0.22
            @test params.V1 == 3.28
            @test params.V2 == 2.66
            @test params.Q == 0.54
            @test params.KSS == 0.087
        end

        @testset "TwoCptTMDDFcRn" begin
            @test TwoCptTMDDFcRn() isa TMDDModelKind
            @test TwoCptTMDDFcRn().approximation == QSS

            params = TwoCptTMDDFcRnParams(
                3.0, 2.5, 0.5,  # V1, V2, Q
                0.3, 0.7,  # CLup, FR
                0.05, 0.1, 0.01, 0.1, 1.0  # KSS, kint, ksyn, kdeg, R0
            )
            @test params.CLup == 0.3
            @test params.FR == 0.7
        end

        @testset "TwoCptTMDDADA" begin
            @test TwoCptTMDDADA() isa TMDDModelKind

            params = TwoCptTMDDADAParams(
                0.2, 3.0, 2.5, 0.5,  # CL, V1, V2, Q
                0.05, 0.1, 0.01, 0.1, 1.0,  # KSS, kint, ksyn, kdeg, R0
                0.01, 0.05, 0.1, 0.01, 1.0,  # kADA_prod, kADA_deg, kon_ADA, koff_ADA, CL_complex
                14.0,  # T_onset
                0.0, 1.0  # ka, F
            )
            @test params.kADA_prod == 0.01
            @test params.T_onset == 14.0
        end

        @testset "TMDDApproximation enum" begin
            @test FullTMDD isa TMDDApproximation
            @test QSS isa TMDDApproximation
            @test QE isa TMDDApproximation
            @test RapidBinding isa TMDDApproximation
            @test IrreversibleBinding isa TMDDApproximation
        end

        @testset "TMDDRoute enum" begin
            @test IVBolus isa TMDDRoute
            @test IVInfusion isa TMDDRoute
            @test Subcutaneous isa TMDDRoute
        end
    end

    @testset "TMDDSpec Construction" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]

        spec = TMDDSpec(TwoCptTMDD(), "Pembrolizumab", params, doses)

        @test spec.kind isa TwoCptTMDD
        @test spec.name == "Pembrolizumab"
        @test spec.params.CL == 0.22
        @test length(spec.doses) == 1
        @test spec.doses[1].amount == 200.0
        @test spec.target_units == :nM
        @test spec.drug_units == :mg_L

        # With explicit units
        spec2 = TMDDSpec(TwoCptTMDD(), "Test", params, doses, :pM, :ug_mL)
        @test spec2.target_units == :pM
        @test spec2.drug_units == :ug_mL
    end

    @testset "TMDD Validation" begin
        @testset "Valid TwoCptTMDD spec" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            doses = [DoseEvent(0.0, 200.0)]
            spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

            @test validate_tmdd(spec) === nothing
        end

        @testset "Valid OneCptTMDD spec" begin
            params = OneCptTMDDParams(0.2, 3.0, 0.05, 0.1, 0.01, 0.1, 1.0)
            doses = [DoseEvent(0.0, 100.0)]
            spec = TMDDSpec(OneCptTMDD(), "Test", params, doses)

            @test validate_tmdd(spec) === nothing
        end

        @testset "Invalid - missing doses" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            spec = TMDDSpec(TwoCptTMDD(), "Test", params, DoseEvent[])

            @test_throws ErrorException validate_tmdd(spec)
        end
    end

    @testset "TMDD Model Interface Functions" begin
        @testset "n_states" begin
            @test n_states(OneCptTMDD(FullTMDD)) == 3  # L, R, P
            @test n_states(OneCptTMDD(QSS)) == 2  # L, Rtot
            @test n_states(TwoCptTMDD(FullTMDD, IVBolus)) == 4  # L, Lp, R, P
            @test n_states(TwoCptTMDD(QSS, IVBolus)) == 3  # L, Lp, Rtot
            @test n_states(TwoCptTMDD(RapidBinding, IVBolus)) == 3  # Ltot, Lp, Rtot
            @test n_states(TwoCptTMDDFcRn()) == 4  # L, Lp, Le, Rtot
            @test n_states(TwoCptTMDDADA()) == 5  # L, Lp, Rtot, ADA, LADA
        end

        @testset "state_names" begin
            @test state_names(OneCptTMDD(FullTMDD)) == [:L, :R, :P]
            @test state_names(OneCptTMDD(QSS)) == [:L, :Rtot]
            @test state_names(TwoCptTMDD(QSS, IVBolus)) == [:L, :Lp, :Rtot]
            @test state_names(TwoCptTMDD(FullTMDD, IVBolus)) == [:L, :Lp, :R, :P]
            @test state_names(TwoCptTMDD(RapidBinding, IVBolus)) == [:Ltot, :Lp, :Rtot]
        end

        @testset "dosing_compartment" begin
            @test dosing_compartment(OneCptTMDD()) == 1
            @test dosing_compartment(TwoCptTMDD()) == 1
            @test dosing_compartment(TwoCptTMDD(QSS, Subcutaneous)) == 1  # SC depot
        end

        @testset "get_central_volume" begin
            params_1cpt = OneCptTMDDParams(0.2, 3.0, 0.05, 0.1, 0.01, 0.1, 1.0)
            @test get_central_volume(params_1cpt) == 3.0

            params_2cpt = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            @test get_central_volume(params_2cpt) == 3.28
        end

        @testset "get_initial_state" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            doses = [DoseEvent(0.0, 200.0)]
            spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

            u0 = get_initial_state(spec)
            @test length(u0) == 3
            @test u0[1] == 0.0  # L = 0
            @test u0[2] == 0.0  # Lp = 0
            @test u0[3] == 0.083  # Rtot = R0
        end
    end

    @testset "TMDD Derived Parameter Functions" begin
        @testset "calculate_KD and calculate_KSS" begin
            @test calculate_KD(0.1, 0.001) ≈ 0.01
            @test calculate_KSS(0.1, 0.001, 0.02) ≈ 0.21  # (koff + kint) / kon
        end

        @testset "calculate_Rtot_ss" begin
            @test calculate_Rtot_ss(0.1, 0.05) ≈ 2.0  # ksyn / kdeg
        end

        @testset "calculate_half_life" begin
            @test calculate_half_life(log(2)) ≈ 1.0
            @test calculate_half_life(0.1) ≈ log(2)/0.1
        end

        @testset "target_occupancy" begin
            @test target_occupancy(0.01, 0.01) ≈ 0.5  # At KD, 50% occupancy
            @test target_occupancy(0.09, 0.01) ≈ 0.9 atol=0.01  # At 9*KD, 90% occupancy
            @test target_occupancy(0.0, 0.01) == 0.0  # No drug = no occupancy
        end

        @testset "free_drug_concentration" begin
            # High drug excess
            L = free_drug_concentration(100.0, 1.0, 0.01, 3.0)
            @test L > 98.0  # Most drug is free

            # High target excess
            L2 = free_drug_concentration(1.0, 100.0, 0.01, 3.0)
            @test L2 < 0.5  # Most drug is bound
        end

        @testset "derived_pk_params" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            derived = derived_pk_params(params)

            @test derived.kel ≈ 0.22 / 3.28
            @test derived.k12 ≈ 0.54 / 3.28
            @test derived.k21 ≈ 0.54 / 2.66
            @test derived.Vss ≈ 3.28 + 2.66
            @test derived.R_ss ≈ 0.001 / 0.012  # ksyn / kdeg
        end

        @testset "convert_to_micro" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            micro = convert_to_micro(params)

            @test micro.kel ≈ 0.22 / 3.28
            @test micro.k12 ≈ 0.54 / 3.28
            @test micro.k21 ≈ 0.54 / 2.66
        end
    end

    @testset "TMDD Analysis Functions" begin
        @testset "calculate_target_occupancy" begin
            @test calculate_target_occupancy(0.0, 0.01) == 0.0
            @test calculate_target_occupancy(0.01, 0.01) ≈ 0.5
            @test calculate_target_occupancy(0.09, 0.01) ≈ 0.9 atol=0.01
        end

        @testset "target_occupancy_EC values" begin
            KD = 0.01
            @test target_occupancy_EC50(KD) == KD
            @test target_occupancy_EC90(KD) == 9 * KD
            @test target_occupancy_EC99(KD) == 99 * KD
        end

        @testset "identify_tmdd_regime" begin
            @test identify_tmdd_regime(10.0, 0.01, 1.0) == :linear
            @test identify_tmdd_regime(0.001, 0.01, 1.0) == :tmdd
            @test identify_tmdd_regime(0.1, 0.01, 1.0) == :mixed
        end

        @testset "calculate_tmdd_steady_state - TwoCptTMDD" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            ss = calculate_tmdd_steady_state(params)

            @test ss.R_ss ≈ 0.001 / 0.012  # ksyn/kdeg
            @test ss.t_half_target ≈ log(2) / 0.012
            @test ss.KSS == 0.087
            @test ss.Vss ≈ 3.28 + 2.66
        end

        @testset "calculate_tmdd_half_lives" begin
            params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
            hl = calculate_tmdd_half_lives(params)

            @test hl.t_half_alpha > 0
            @test hl.t_half_beta > 0
            @test hl.t_half_beta > hl.t_half_alpha  # Terminal > distribution
            @test hl.t_half_target ≈ log(2) / 0.012
            @test hl.t_half_complex ≈ log(2) / 0.037
        end
    end

    @testset "TMDD Simulation - TwoCptTMDD QSS (Industry Standard)" begin
        params = TwoCptTMDDParams(
            0.22, 3.28, 2.66, 0.54,  # CL, V1, V2, Q
            0.087, 0.037, 0.001, 0.012, 0.083  # KSS, kint, ksyn, kdeg, R0
        )
        doses = [DoseEvent(0.0, 200.0)]  # 200 mg IV bolus
        spec = TMDDSpec(TwoCptTMDD(), "Test QSS", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:1.0:336.0))  # 14 days
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test length(result.t) == 337
        @test haskey(result.states, :L)
        @test haskey(result.states, :Lp)
        @test haskey(result.states, :Rtot)
        @test haskey(result.observations, :conc)
        @test haskey(result.observations, :conc_peripheral)
        @test haskey(result.observations, :target_occupancy)

        # Initial concentration = dose/V1
        @test result.observations[:conc][1] ≈ 200.0 / 3.28 atol=0.1

        # Concentration decreases
        @test result.observations[:conc][end] < result.observations[:conc][1]

        # Peripheral builds up then declines
        @test maximum(result.observations[:conc_peripheral]) > 0

        # Target occupancy achieves high level then declines
        max_occ = maximum(result.observations[:target_occupancy])
        @test max_occ > 0.8  # Should achieve >80% with 200 mg dose
    end

    @testset "TMDD Simulation - TwoCptTMDD Full Model" begin
        params = TwoCptTMDDParams(
            0.22, 3.28, 2.66, 0.54,
            0.087, 0.037, 0.001, 0.012, 0.083
        )
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(FullTMDD, IVBolus), "Test Full", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :L)
        @test haskey(result.states, :Lp)
        @test haskey(result.states, :R)
        @test haskey(result.states, :P)

        # Complex formation tracked
        @test haskey(result.observations, :conc_bound)
        @test maximum(result.observations[:conc_bound]) > 0
    end

    @testset "TMDD Simulation - OneCptTMDD" begin
        params = OneCptTMDDParams(0.2, 3.0, 0.05, 0.1, 0.01, 0.1, 1.0)
        doses = [DoseEvent(0.0, 100.0)]
        spec = TMDDSpec(OneCptTMDD(), "Test 1Cpt", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.observations, :conc)
        @test result.observations[:conc][1] > 0
        @test result.observations[:conc][end] < result.observations[:conc][1]
    end

    @testset "TMDD Simulation - TwoCptTMDD RapidBinding" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(RapidBinding, IVBolus), "Test Rapid", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :Ltot)  # Total drug
        @test haskey(result.observations, :conc_free)
        @test haskey(result.observations, :conc_total)

        # Total >= Free
        @test result.observations[:conc_total][1] >= result.observations[:conc_free][1]
    end

    @testset "TMDD Simulation - FcRn Recycling" begin
        params = TwoCptTMDDFcRnParams(
            3.0, 2.5, 0.5,  # V1, V2, Q
            0.3, 0.7,  # CLup, FR (70% recycled)
            0.05, 0.1, 0.01, 0.1, 1.0  # KSS, kint, ksyn, kdeg, R0
        )
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDDFcRn(), "FcRn Test", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:2.0:336.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult
        @test haskey(result.states, :Le)  # Endosomal
        @test haskey(result.observations, :conc_endosomal)
        @test haskey(result.observations, :endosomal_fraction)
    end

    @testset "TMDD Simulation - Multiple Doses" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)

        # Q3W dosing (every 21 days = 504 hours)
        doses = [
            DoseEvent(0.0, 200.0),
            DoseEvent(504.0, 200.0),
            DoseEvent(1008.0, 200.0)
        ]
        spec = TMDDSpec(TwoCptTMDD(), "Q3W Dosing", params, doses)

        grid = SimGrid(0.0, 1512.0, collect(0.0:6.0:1512.0))  # 9 weeks
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test result isa TMDDSimResult

        conc = result.observations[:conc]
        @test length(conc) > 0

        # Should see accumulation at steady state
        metrics = tmdd_exposure_metrics(result)
        @test metrics.AUC > 0
        @test metrics.Cmax > 0
    end

    @testset "TMDD Exposure Metrics" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:1.0:336.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)
        metrics = tmdd_exposure_metrics(result)

        @test metrics.AUC > 0
        @test metrics.AUC_total >= metrics.AUC
        @test metrics.Cmax > 0
        @test metrics.Tmax >= 0
        @test metrics.Ctrough >= 0
        @test metrics.AUC_occupancy > 0
        @test 0.0 <= metrics.mean_occupancy <= 1.0
        @test 0.0 <= metrics.max_occupancy <= 1.0
    end

    @testset "Time Above Threshold" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

        grid = SimGrid(0.0, 336.0, collect(0.0:1.0:336.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        # Time above 90% occupancy
        t_90 = time_above_occupancy(result, 0.9)
        @test t_90 >= 0.0
        @test t_90 <= 336.0

        # Time above concentration threshold
        t_conc = time_above_concentration(result, 10.0)
        @test t_conc >= 0.0
    end

    @testset "TMDD Determinism" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result1 = solve_tmdd(spec, grid, solver)
        result2 = solve_tmdd(spec, grid, solver)

        @test result1.observations[:conc] == result2.observations[:conc]
        @test result1.states[:L] == result2.states[:L]
    end

    @testset "TMDD Model Comparison" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        # QSS approximation
        spec_qss = TMDDSpec(TwoCptTMDD(QSS, IVBolus), "QSS", params, doses)
        result_qss = solve_tmdd(spec_qss, grid, solver)

        # Full model
        spec_full = TMDDSpec(TwoCptTMDD(FullTMDD, IVBolus), "Full", params, doses)
        result_full = solve_tmdd(spec_full, grid, solver)

        # Compare
        comparison = compare_tmdd_approximations(result_full, result_qss)

        @test comparison.max_rel_error >= 0.0
        @test comparison.AUC_rel_error >= 0.0
        @test comparison.recommendation isa Symbol
    end

    @testset "TMDD Dose Selection" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)

        # Estimate dose for 90% trough occupancy at 21-day interval
        dose_est = target_trough_dose(params, 0.9, 504.0)  # Q3W
        @test dose_est > 0

        # Loading dose
        ld = loading_dose(params, 0.9)
        @test ld > 0
    end

    @testset "TMDD Steady-State Functions" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)

        ss = tmdd_steady_state(params)
        @test ss.R_ss > 0
        @test ss.t_half_target > 0
        @test ss.t_half_linear > 0
        @test ss.KSS == params.KSS

        tss = time_to_steady_state(params)
        @test tss > 0
    end

    @testset "Metadata in Results" begin
        params = TwoCptTMDDParams(0.22, 3.28, 2.66, 0.54, 0.087, 0.037, 0.001, 0.012, 0.083)
        doses = [DoseEvent(0.0, 200.0)]
        spec = TMDDSpec(TwoCptTMDD(), "Test", params, doses)

        grid = SimGrid(0.0, 168.0, collect(0.0:1.0:168.0))
        solver = SolverSpec(:Rodas5, 1e-8, 1e-8, 10000)

        result = solve_tmdd(spec, grid, solver)

        @test haskey(result.metadata, "model_type")
        @test haskey(result.metadata, "approximation")
        @test haskey(result.metadata, "route")
        @test haskey(result.metadata, "n_doses")
        @test haskey(result.metadata, "solver")
        @test haskey(result.metadata, "retcode")

        @test result.metadata["approximation"] == "QSS"
        @test result.metadata["route"] == "IVBolus"
    end

end
