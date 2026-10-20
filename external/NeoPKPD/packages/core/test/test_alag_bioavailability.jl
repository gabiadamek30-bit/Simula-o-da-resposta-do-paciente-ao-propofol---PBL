# Test suite for ALAG (Absorption Lag Time) and Bioavailability (F) implementation
# Industry-standard tests based on NONMEM/Monolix validation scenarios

using Test
using NeoPKPD

@testset "ALAG and Bioavailability Support" begin

    @testset "AbsorptionModifiers construction and validation" begin
        # Default modifiers
        mod_default = AbsorptionModifiers()
        @test mod_default.alag == 0.0
        @test mod_default.bioavailability == 1.0

        # Positional constructor
        mod_positional = AbsorptionModifiers(1.0, 0.7)
        @test mod_positional.alag == 1.0
        @test mod_positional.bioavailability == 0.7
    end

    @testset "ModelSpecWithModifiers construction" begin
        # Create model spec with modifiers
        spec = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Test Model",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),  # Ka, CL, V
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.5, 0.8)  # alag, F
        )

        @test spec.kind isa OneCompOralFirstOrder
        @test spec.name == "Test Model"
        @test spec.params.Ka == 1.5
        @test spec.absorption_modifiers.alag == 0.5
        @test spec.absorption_modifiers.bioavailability == 0.8

        # with_modifiers helper
        base_spec = ModelSpec(
            OneCompOralFirstOrder(),
            "Base Model",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)]
        )

        converted = with_modifiers(base_spec, AbsorptionModifiers(0.3, 0.9))
        @test converted.absorption_modifiers.alag == 0.3
        @test converted.absorption_modifiers.bioavailability == 0.9

        # has_absorption_modifiers checks
        @test has_absorption_modifiers(base_spec) == false
        @test has_absorption_modifiers(spec) == true

        spec_default_mods = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Default Mods",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers()  # default: alag=0, F=1
        )
        @test has_absorption_modifiers(spec_default_mods) == false
    end

    @testset "OneCompOralFirstOrder with ALAG" begin
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Reference: no ALAG
        spec_no_alag = ModelSpec(
            OneCompOralFirstOrder(),
            "No ALAG",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),  # Ka, CL, V
            [DoseEvent(0.0, 300.0)]
        )
        result_no_alag = simulate(spec_no_alag, grid, solver)

        # With ALAG = 0.5 hours
        spec_with_alag = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "With ALAG",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.5, 1.0)  # ALAG=0.5, F=1.0
        )
        result_with_alag = simulate(spec_with_alag, grid, solver)

        # Verify ALAG shifts the profile
        # At t=0, with ALAG the concentration should be 0 (dose hasn't started absorbing)
        @test result_with_alag.observations[:conc][1] ≈ 0.0 atol=1e-10

        # The profile should be shifted by 0.5 hours
        idx_05 = findfirst(t -> t ≈ 0.5, result_with_alag.t)
        @test result_with_alag.observations[:conc][idx_05] ≈ 0.0 atol=1e-6  # Dose just started

        # At t=1.0 with ALAG, should match t=0.5 without ALAG (approximately)
        idx_1 = findfirst(t -> t ≈ 1.0, result_with_alag.t)
        idx_05_no_alag = findfirst(t -> t ≈ 0.5, result_no_alag.t)
        @test result_with_alag.observations[:conc][idx_1] ≈ result_no_alag.observations[:conc][idx_05_no_alag] rtol=0.05

        # Check metadata contains ALAG info
        @test haskey(result_with_alag.metadata, "absorption_modifiers")
        @test result_with_alag.metadata["absorption_modifiers"]["ALAG"] == 0.5
    end

    @testset "OneCompOralFirstOrder with Bioavailability" begin
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Reference: F=1.0 (100% bioavailability)
        spec_f100 = ModelSpec(
            OneCompOralFirstOrder(),
            "F=100%",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)]
        )
        result_f100 = simulate(spec_f100, grid, solver)

        # With F=0.5 (50% bioavailability)
        spec_f50 = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "F=50%",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.0, 0.5)  # ALAG=0, F=0.5
        )
        result_f50 = simulate(spec_f50, grid, solver)

        # Concentrations should be exactly half with F=0.5
        for i in 1:length(result_f100.t)
            @test result_f50.observations[:conc][i] ≈ result_f100.observations[:conc][i] * 0.5 rtol=1e-6
        end

        # Check metadata
        @test result_f50.metadata["absorption_modifiers"]["F"] == 0.5
    end

    @testset "Combined ALAG and Bioavailability" begin
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Combined: ALAG=0.5, F=0.7
        spec_combined = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Combined",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.5, 0.7)
        )
        result_combined = simulate(spec_combined, grid, solver)

        # At t=0, concentration should be 0 (due to ALAG)
        @test result_combined.observations[:conc][1] ≈ 0.0 atol=1e-10

        # AUC should be proportional to F (for linear PK)
        # Simple check: Cmax should be ~70% of F=100% case
        spec_f100 = ModelSpec(
            OneCompOralFirstOrder(),
            "F=100%",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)]
        )
        result_f100 = simulate(spec_f100, grid, solver)

        cmax_f100 = maximum(result_f100.observations[:conc])
        cmax_combined = maximum(result_combined.observations[:conc])
        @test cmax_combined ≈ cmax_f100 * 0.7 rtol=0.02  # Allow 2% tolerance for lag time effect
    end

    @testset "TwoCompOral with ALAG and F" begin
        grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Two-compartment oral with ALAG and F
        spec = ModelSpecWithModifiers(
            TwoCompOral(),
            "TwoComp ALAG+F",
            TwoCompOralParams(2.0, 3.0, 20.0, 5.0, 40.0),  # Ka, CL, V1, Q, V2
            [DoseEvent(0.0, 500.0)],
            AbsorptionModifiers(1.0, 0.65)
        )
        result = simulate(spec, grid, solver)

        # Verify zero concentration at t<ALAG
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10  # t=0

        # Check all compartments are present
        @test haskey(result.states, :A_gut)
        @test haskey(result.states, :A_central)
        @test haskey(result.states, :A_peripheral)

        # Verify metadata
        @test result.metadata["absorption_modifiers"]["ALAG"] == 1.0
        @test result.metadata["absorption_modifiers"]["F"] == 0.65
    end

    @testset "TransitAbsorption with ALAG and F" begin
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Transit absorption with ALAG (stacked delays)
        spec = ModelSpecWithModifiers(
            TransitAbsorption(),
            "Transit ALAG+F",
            TransitAbsorptionParams(3, 2.0, 1.5, 2.0, 30.0),  # N, Ktr, Ka, CL, V
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.25, 0.85)
        )
        result = simulate(spec, grid, solver)

        # Verify ALAG delays the transit chain start
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10

        # Check transit compartments exist
        @test haskey(result.states, :Transit_1)
        @test haskey(result.states, :Transit_2)
        @test haskey(result.states, :Transit_3)

        @test result.metadata["N_transit"] == 3
    end

    @testset "IV models with ALAG (delayed infusion start)" begin
        grid = SimGrid(0.0, 12.0, collect(0.0:0.25:12.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # OneCompIVBolus with ALAG (delayed IV, e.g., start of infusion)
        spec = ModelSpecWithModifiers(
            OneCompIVBolus(),
            "Delayed IV",
            OneCompIVBolusParams(2.0, 30.0),  # CL, V
            [DoseEvent(0.0, 100.0)],
            AbsorptionModifiers(0.5, 1.0)  # 30 min delay
        )
        result = simulate(spec, grid, solver)

        # At t=0, concentration should be 0 (IV hasn't started)
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10

        # At t=0.5+ε, concentration should jump
        idx_post_alag = findfirst(t -> t > 0.5, result.t)
        @test result.observations[:conc][idx_post_alag] > 0.0
    end

    @testset "Multiple dose ALAG/F handling" begin
        grid = SimGrid(0.0, 48.0, collect(0.0:1.0:48.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Multiple doses with ALAG
        spec = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Multiple Doses ALAG",
            OneCompOralFirstOrderParams(1.0, 2.0, 30.0),
            [
                DoseEvent(0.0, 200.0),
                DoseEvent(12.0, 200.0),
                DoseEvent(24.0, 200.0)
            ],
            AbsorptionModifiers(0.5, 0.8)
        )
        result = simulate(spec, grid, solver)

        # Each dose should be delayed by ALAG
        # At t=0, should still be 0
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10

        # Concentration should be non-zero after first ALAG
        idx_1 = findfirst(t -> t >= 1.0, result.t)
        @test result.observations[:conc][idx_1] > 0.0

        # Check concentrations rise after each dose+ALAG
        idx_13 = findfirst(t -> t >= 13.0, result.t)
        idx_12 = findfirst(t -> t >= 12.0, result.t)
        @test result.observations[:conc][idx_13] > 0.0
    end

    @testset "NONMEM-style validation: Theophylline oral" begin
        # Based on NONMEM ADVAN2 theophylline example
        # Typical values: Ka=1.5/hr, CL=2.0 L/hr, V=30 L, ALAG=0.4 hr, F=0.9
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Standard oral theophylline dosing
        spec = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Theophylline",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],  # 300 mg oral dose
            AbsorptionModifiers(0.4, 0.9)  # ALAG=0.4, F=0.9
        )
        result = simulate(spec, grid, solver)

        # Expected characteristics:
        # 1. Zero concentration until t ≈ ALAG
        # 2. Tmax should be shifted by ALAG
        # 3. AUC should be ~90% of IV equivalent

        # Check zero before ALAG
        @test result.observations[:conc][1] ≈ 0.0 atol=1e-10

        # Find Cmax and Tmax
        conc = result.observations[:conc]
        cmax_idx = argmax(conc)
        tmax = result.t[cmax_idx]
        cmax = conc[cmax_idx]

        # Tmax should be > ALAG
        @test tmax > 0.4

        # For linear 1-comp, analytical Cmax with ALAG is at:
        # tmax = ALAG + log(Ka/kel) / (Ka - kel) where kel = CL/V
        kel = 2.0 / 30.0
        ka = 1.5
        analytical_tmax = 0.4 + log(ka / kel) / (ka - kel)
        @test tmax ≈ analytical_tmax rtol=0.1  # 10% tolerance

        # Verify Cmax is reasonable (should be < dose/V due to F<1 and distribution)
        @test cmax < 300.0 * 0.9 / 30.0  # Must be less than F*Dose/V
        @test cmax > 0.0  # Must be positive
    end

    @testset "Edge cases" begin
        grid = SimGrid(0.0, 12.0, collect(0.0:0.5:12.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # ALAG longer than first observation interval
        spec_long_alag = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Long ALAG",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(3.0, 1.0)  # 3 hour lag
        )
        result_long = simulate(spec_long_alag, grid, solver)

        # Should be zero for t < 3
        for i in 1:length(result_long.t)
            if result_long.t[i] < 3.0
                @test result_long.observations[:conc][i] ≈ 0.0 atol=1e-10
            end
        end

        # Very low bioavailability
        spec_low_f = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Low F",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)],
            AbsorptionModifiers(0.0, 0.05)  # 5% bioavailability
        )
        result_low_f = simulate(spec_low_f, grid, solver)

        # Reference without modifiers
        spec_ref = ModelSpec(
            OneCompOralFirstOrder(),
            "Reference",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(0.0, 300.0)]
        )
        result_ref = simulate(spec_ref, grid, solver)

        # Concentrations should be 5% of reference
        for i in 1:length(result_low_f.t)
            @test result_low_f.observations[:conc][i] ≈ result_ref.observations[:conc][i] * 0.05 rtol=1e-5
        end
    end

    @testset "Dose at non-zero time with ALAG" begin
        grid = SimGrid(0.0, 24.0, collect(0.0:0.5:24.0))
        solver = SolverSpec(:Tsit5, 1e-8, 1e-8, 10000)

        # Dose at t=6 with ALAG=1.0 -> effective dose time = 7.0
        spec = ModelSpecWithModifiers(
            OneCompOralFirstOrder(),
            "Delayed Start",
            OneCompOralFirstOrderParams(1.5, 2.0, 30.0),
            [DoseEvent(6.0, 300.0)],  # Dose at t=6
            AbsorptionModifiers(1.0, 1.0)
        )
        result = simulate(spec, grid, solver)

        # Should be zero until t=7 (dose_time + ALAG)
        for i in 1:length(result.t)
            if result.t[i] < 7.0
                @test result.observations[:conc][i] ≈ 0.0 atol=1e-10
            end
        end

        # Should have concentration after t=7
        idx_8 = findfirst(t -> t >= 8.0, result.t)
        @test result.observations[:conc][idx_8] > 0.0
    end
end
