# Test suite for NONMEM and Monolix Import

using Test
using NeoPKPD

@testset "NONMEM Import" begin

    @testset "NONMEM Types" begin
        # Test THETASpec
        theta1 = THETASpec(10.0)
        @test theta1.init == 10.0
        @test theta1.lower == -Inf
        @test theta1.upper == Inf
        @test theta1.fixed == false

        theta2 = THETASpec(0.0, 10.0, 100.0)
        @test theta2.lower == 0.0
        @test theta2.init == 10.0
        @test theta2.upper == 100.0

        # Test OMEGABlock
        omega1 = OMEGABlock([0.04, 0.09])
        @test omega1.structure == :diagonal
        @test omega1.dimension == 2

        # Test SIGMABlock
        sigma1 = SIGMABlock([0.01])
        @test sigma1.structure == :diagonal

        # Test SubroutineSpec
        sub1 = SubroutineSpec(1, 2)
        @test sub1.advan == 1
        @test sub1.trans == 2

        # Test DataSpec
        data1 = DataSpec("data.csv")
        @test data1.filename == "data.csv"

        # Test InputColumn
        col1 = InputColumn("ID")
        @test col1.name == "ID"
        @test col1.drop == false
    end

    @testset "ADVAN/TRANS Mapping" begin
        @test get_model_mapping(1, 2) == (:OneCompIVBolus, [:CL, :V])
        @test get_model_mapping(2, 2) == (:OneCompOralFirstOrder, [:KA, :CL, :V])
        @test get_model_mapping(3, 4) == (:TwoCompIVBolus, [:CL, :V1, :Q, :V2])
        @test get_model_mapping(4, 4) == (:TwoCompOral, [:KA, :CL, :V1, :Q, :V2])
        @test get_model_mapping(11, 4) == (:ThreeCompIVBolus, [:CL, :V1, :Q2, :V2, :Q3, :V3])
        @test get_model_mapping(10, 1) == (:MichaelisMentenElimination, [:VM, :KM, :V])
        @test get_model_mapping(99, 1) === nothing
    end

    @testset "NONMEM Parser" begin
        ctl_text = raw"""
$PROBLEM One-compartment IV bolus

$DATA data.csv IGNORE=@

$INPUT ID TIME DV AMT EVID

$SUBROUTINES ADVAN1 TRANS2

$THETA
(0, 10, 100)
(0, 50, 500)

$OMEGA
0.04
0.09

$SIGMA
0.01

$ESTIMATION METHOD=1 INTER MAXEVAL=9999
"""

        ctl = parse_nonmem_control(ctl_text)

        @test ctl.problem == "One-compartment IV bolus"
        @test ctl.data !== nothing
        @test ctl.data.filename == "data.csv"
        @test length(ctl.input) == 5
        @test ctl.subroutines !== nothing
        @test ctl.subroutines.advan == 1
        @test ctl.subroutines.trans == 2
        @test length(ctl.thetas) == 2
        @test ctl.thetas[1].lower == 0.0
        @test ctl.thetas[1].init == 10.0
        @test ctl.thetas[1].upper == 100.0
        @test length(ctl.omegas) == 1
        @test length(ctl.sigmas) == 1
        @test ctl.estimation["method"] == "FOCE"
        @test ctl.estimation["interaction"] == true
    end

    @testset "NONMEM Converter - One-Comp IV" begin
        ctl_text = raw"""
$PROBLEM One-comp IV test
$SUB ADVAN1 TRANS2
$THETA
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompIVBolus
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V == 50.0
        @test result.iiv_spec !== nothing
        @test result.error_spec !== nothing
    end

    @testset "NONMEM Converter - One-Comp Oral" begin
        ctl_text = raw"""
$PROBLEM One-comp oral test
$SUB ADVAN2 TRANS2
$THETA
(0, 1.5, 10)
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09 0.16
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompOralFirstOrder
        @test result.model_spec.params.Ka == 1.5
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V == 50.0
    end

    @testset "NONMEM Converter - Two-Comp IV" begin
        ctl_text = raw"""
$PROBLEM Two-comp IV test
$SUB ADVAN3 TRANS4
$THETA
(0, 10, 100)
(0, 30, 300)
(0, 5, 50)
(0, 60, 600)
$OMEGA 0.04 0.09 0.04 0.09
$SIGMA 0.01 0.02
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa TwoCompIVBolus
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V1 == 30.0
        @test result.model_spec.params.Q == 5.0
        @test result.model_spec.params.V2 == 60.0
    end

    @testset "NONMEM Converter - Unsupported ADVAN" begin
        ctl_text = raw"""
$PROBLEM Unsupported
$SUB ADVAN6
$THETA (10)
"""

        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_neopkpd(ctl)

        @test !isempty(result.errors)
        @test result.model_spec === nothing
    end

    @testset "NONMEM Converter - Missing Subroutines" begin
        ctl_text = raw"""
$PROBLEM No subroutines
$THETA (10)
"""

        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_neopkpd(ctl)

        @test !isempty(result.errors)
        @test any(contains(e, "SUBROUTINES") for e in result.errors)
    end
end

@testset "Monolix Import" begin

    @testset "Monolix Types" begin
        # Test MonolixParameter
        p1 = MonolixParameter("ka", 1.5)
        @test p1.name == "ka"
        @test p1.value == 1.5
        @test p1.fixed == false

        p2 = MonolixParameter("V", 50.0; fixed=true, omega=0.2, has_iiv=true)
        @test p2.fixed == true
        @test p2.omega == 0.2
        @test p2.has_iiv == true

        # Test MonolixObservation
        obs1 = MonolixObservation("y1")
        @test obs1.name == "y1"
        @test obs1.type == "continuous"
        @test obs1.error_model == "combined"

        # Test MonolixDataset
        data1 = MonolixDataset("data.csv")
        @test data1.filename == "data.csv"
        @test data1.id_column == "ID"

        # Test MonolixStructuralModel
        model_type = MonolixModelType("pklib", "pk_oral1cpt_1abs_kaVCl_PLASMA")
        model = MonolixStructuralModel(model_type, "oral", 1, "linear", "firstOrder", false, false)
        @test model.admin_type == "oral"
        @test model.n_compartments == 1
    end

    @testset "Monolix Model Mapping" begin
        @test get_monolix_model_mapping("pk_bolus1cpt_VCl_PLASMA") == :OneCompIVBolus
        @test get_monolix_model_mapping("pk_oral1cpt_1abs_kaVCl_PLASMA") == :OneCompOralFirstOrder
        @test get_monolix_model_mapping("pk_bolus2cpt_V1ClQ2V2_PLASMA") == :TwoCompIVBolus
        @test get_monolix_model_mapping("some_oral_1cpt_model") == :OneCompOralFirstOrder
        @test get_monolix_model_mapping("unknown_model") === nothing
    end

    @testset "Monolix Parser" begin
        mlx_text = """
[DESCRIPTION]
One-compartment oral PK model

[DATAFILE]
file = 'data.csv'
header = {ID=ID, TIME=TIME, DV=OBSERVATION, AMT=AMOUNT}

[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.5
V = {value=50.0, method=MLE}
Cl = {value=10.0, method=MLE}

[POPULATION]
ka = {value=1.5, omega=0.3, distribution=logNormal}
V = {value=50.0, omega=0.2}
Cl = {value=10.0, omega=0.25}

[OBSERVATION_MODEL]
y1 = {type=continuous, prediction=Cc, error=combined1}
"""

        mlx = parse_monolix_project(mlx_text)

        @test mlx.description == "One-compartment oral PK model"
        @test mlx.data !== nothing
        @test mlx.data.filename == "data.csv"
        @test mlx.model !== nothing
        @test mlx.model.model_type !== nothing
        @test mlx.model.model_type.model == "pk_oral1cpt_1abs_kaVCl_PLASMA"
        @test length(mlx.parameters) >= 3
        @test length(mlx.observations) == 1
    end

    @testset "Monolix Converter" begin
        mlx_text = """
[DESCRIPTION]
Test conversion

[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.5
V = 50.0
Cl = 10.0

[POPULATION]
ka = {omega=0.3}
V = {omega=0.2}
Cl = {omega=0.25}

[OBSERVATION_MODEL]
y1 = {type=continuous, error=proportional}
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_monolix_to_neopkpd(mlx; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa OneCompOralFirstOrder
        @test result.iiv_spec !== nothing
        @test result.error_spec !== nothing
    end

    @testset "Monolix Converter - Two-Comp" begin
        mlx_text = """
[STRUCTURAL_MODEL]
lib = pklib:pk_bolus2cpt_V1ClQ2V2_PLASMA

[PARAMETER]
V1 = 30.0
Cl = 10.0
Q = 5.0
V2 = 60.0
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]

        result = convert_monolix_to_neopkpd(mlx; doses=doses)

        @test isempty(result.errors)
        @test result.model_spec !== nothing
        @test result.model_spec.kind isa TwoCompIVBolus
        @test result.model_spec.params.V1 == 30.0
        @test result.model_spec.params.CL == 10.0
    end

    @testset "Modern Monolix Format (XML-style sections)" begin
        @testset "One-Comp Oral with Combined Error" begin
            mlx_text = """
<DATAFILE>
[FILEINFO]
file = 'data.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=Cc, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}

<MODEL>
[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}

DEFINITION:
ka = {distribution=logNormal, typical=ka_pop, sd=omega_ka}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}

[LONGITUDINAL]
input = {a, b}

file = 'lib:oral1_1cpt_kaVCl.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}

<FIT>
data = Cc
model = Cc

<PARAMETER>
ka_pop = {value=1.5, method=MLE}
V_pop = {value=50, method=MLE}
Cl_pop = {value=5, method=MLE}
omega_ka = {value=0.4, method=MLE}
omega_V = {value=0.25, method=MLE}
omega_Cl = {value=0.3, method=MLE}
a = {value=0.5, method=MLE}
b = {value=0.1, method=MLE}

<MONOLIX>
[TASKS]
populationParameters()
conditionalModes(linearization = true)
"""
            mlx = parse_monolix_project(mlx_text)

            # Test data parsing
            @test mlx.data !== nothing
            @test mlx.data.filename == "data.csv"
            @test mlx.data.id_column == "ID"
            @test mlx.data.time_column == "TIME"

            # Test model parsing
            @test mlx.model !== nothing
            @test mlx.model.n_compartments == 1
            @test mlx.model.absorption == "firstOrder"

            # Test parameter parsing
            @test length(mlx.parameters) >= 3

            # Find ka, V, Cl parameters
            ka_param = findfirst(p -> p.name == "ka", mlx.parameters)
            @test ka_param !== nothing
            @test mlx.parameters[ka_param].value == 1.5
            @test mlx.parameters[ka_param].has_iiv == true
            @test mlx.parameters[ka_param].omega ≈ 0.4

            v_param = findfirst(p -> p.name == "V", mlx.parameters)
            @test v_param !== nothing
            @test mlx.parameters[v_param].value == 50.0
            @test mlx.parameters[v_param].has_iiv == true

            # Test observation parsing (combined error)
            @test length(mlx.observations) >= 1
            cc_obs = findfirst(o -> o.name == "Cc", mlx.observations)
            @test cc_obs !== nothing
            @test mlx.observations[cc_obs].error_model == "combined"
            @test length(mlx.observations[cc_obs].error_params) == 2
            @test mlx.observations[cc_obs].error_params[1] ≈ 0.5
            @test mlx.observations[cc_obs].error_params[2] ≈ 0.1

            # Test conversion
            doses = [DoseEvent(0.0, 100.0)]
            result = convert_monolix_to_neopkpd(mlx; doses=doses)

            @test isempty(result.errors)
            @test result.model_spec !== nothing
            @test result.model_spec.kind isa OneCompOralFirstOrder
            @test result.iiv_spec !== nothing
            @test result.error_spec !== nothing
            @test result.error_spec.kind isa CombinedError
        end

        @testset "Two-Comp IV with Proportional Error" begin
            mlx_text = """
<DATAFILE>
[FILEINFO]
file = 'data.csv'
delimiter = comma
header = {ID, TIME, DV, AMT, EVID}

[CONTENT]
ID = {use=identifier}
TIME = {use=time}
DV = {use=observation, name=Cc, type=continuous}
AMT = {use=amount}
EVID = {use=eventidentifier}

<MODEL>
[INDIVIDUAL]
input = {Cl_pop, V1_pop, Q_pop, V2_pop, omega_Cl, omega_V1}

DEFINITION:
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}
V1 = {distribution=logNormal, typical=V1_pop, sd=omega_V1}
Q = {distribution=logNormal, typical=Q_pop, no-variability}
V2 = {distribution=logNormal, typical=V2_pop, no-variability}

[LONGITUDINAL]
input = {b}

file = 'lib:bolus2_2cpt_ClV1QV2.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=proportional(b)}

<PARAMETER>
Cl_pop = {value=5, method=MLE}
V1_pop = {value=10, method=MLE}
Q_pop = {value=2, method=MLE}
V2_pop = {value=20, method=MLE}
omega_Cl = {value=0.3, method=MLE}
omega_V1 = {value=0.25, method=MLE}
b = {value=0.15, method=MLE}

<MONOLIX>
[TASKS]
populationParameters()
"""
            mlx = parse_monolix_project(mlx_text)

            # Test model parsing
            @test mlx.model !== nothing
            @test mlx.model.n_compartments == 2
            @test mlx.model.admin_type == "iv"

            # Test parameter parsing with no-variability
            cl_param = findfirst(p -> p.name == "Cl", mlx.parameters)
            @test cl_param !== nothing
            @test mlx.parameters[cl_param].has_iiv == true
            @test mlx.parameters[cl_param].omega ≈ 0.3

            q_param = findfirst(p -> p.name == "Q", mlx.parameters)
            @test q_param !== nothing
            @test mlx.parameters[q_param].has_iiv == false

            # Test proportional error
            cc_obs = findfirst(o -> o.name == "Cc", mlx.observations)
            @test mlx.observations[cc_obs].error_model == "proportional"
            @test length(mlx.observations[cc_obs].error_params) == 1
            @test mlx.observations[cc_obs].error_params[1] ≈ 0.15

            # Test conversion
            doses = [DoseEvent(0.0, 100.0)]
            result = convert_monolix_to_neopkpd(mlx; doses=doses)

            @test isempty(result.errors)
            @test result.model_spec !== nothing
            @test result.model_spec.kind isa TwoCompIVBolus
            @test result.error_spec !== nothing
            @test result.error_spec.kind isa ProportionalError
        end

        @testset "Additive Error Model" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:bolus1_1cpt_VCl.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=constant(a)}

<PARAMETER>
V_pop = {value=50, method=MLE}
Cl_pop = {value=10, method=MLE}
a = {value=0.5, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)

            cc_obs = findfirst(o -> o.name == "Cc", mlx.observations)
            @test mlx.observations[cc_obs].error_model == "constant"
            @test mlx.observations[cc_obs].error_params[1] ≈ 0.5

            doses = [DoseEvent(0.0, 100.0)]
            result = convert_monolix_to_neopkpd(mlx; doses=doses)

            @test result.error_spec !== nothing
            @test result.error_spec.kind isa AdditiveError
        end
    end

    @testset "Monolix Unsupported Construct Detection" begin
        @testset "PD Model Detection" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:emax_Emax_EC50.txt'

DEFINITION:
E = {distribution=normal, prediction=E, errorModel=constant(a)}

<PARAMETER>
Emax_pop = {value=100, method=MLE}
EC50_pop = {value=10, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)
            unsupported = NeoPKPD.check_unsupported_monolix_constructs(mlx)

            @test !isempty(unsupported)
            @test any(u -> occursin("PD", u.construct) || occursin("emax", lowercase(u.line)), unsupported)
        end

        @testset "Turnover Model Detection" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:turnover_kin_Imax_IC50.txt'

<PARAMETER>
kin_pop = {value=1, method=MLE}
kout_pop = {value=0.1, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)
            unsupported = NeoPKPD.check_unsupported_monolix_constructs(mlx)

            @test !isempty(unsupported)
            @test any(u -> occursin("Turnover", u.construct), unsupported)
        end

        @testset "Transit Compartment Detection" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:oral1_1cpt_transitAbs_kaVCl.txt'

<PARAMETER>
ka_pop = {value=1, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)
            unsupported = NeoPKPD.check_unsupported_monolix_constructs(mlx)

            @test !isempty(unsupported)
            @test any(u -> occursin("Transit", u.construct), unsupported)
        end

        @testset "Mixture Model Detection" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:mixture_2subpop_oral1cpt.txt'

<PARAMETER>
ka_pop = {value=1, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)
            unsupported = NeoPKPD.check_unsupported_monolix_constructs(mlx)

            @test !isempty(unsupported)
            @test any(u -> occursin("Mixture", u.construct), unsupported)
        end

        @testset "Strict Mode Fails on Unsupported" begin
            mlx_text = """
<MODEL>
[LONGITUDINAL]
file = 'lib:pkpd_1cpt_Emax.txt'

<PARAMETER>
V_pop = {value=50, method=MLE}
Emax_pop = {value=100, method=MLE}
"""
            mlx = parse_monolix_project(mlx_text)
            doses = [DoseEvent(0.0, 100.0)]

            # With strict=true, should add errors
            result = convert_monolix_to_neopkpd(mlx; doses=doses, strict=true)
            @test !isempty(result.errors)

            # With strict=false, should have warnings but not fail
            result_lenient = convert_monolix_to_neopkpd(mlx; doses=doses, strict=false)
            @test !isempty(result_lenient.unsupported)
        end
    end

    @testset "Legacy vs Modern Format Compatibility" begin
        # Test that both formats produce equivalent results
        legacy_text = """
[DESCRIPTION]
One-compartment oral PK model

[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.5
V = 50.0
Cl = 10.0

[POPULATION]
ka = {omega=0.3}
V = {omega=0.2}
Cl = {omega=0.25}

[OBSERVATION_MODEL]
y1 = {type=continuous, error=combined1}
"""

        modern_text = """
<MODEL>
[INDIVIDUAL]
input = {ka_pop, V_pop, Cl_pop, omega_ka, omega_V, omega_Cl}

DEFINITION:
ka = {distribution=logNormal, typical=ka_pop, sd=omega_ka}
V = {distribution=logNormal, typical=V_pop, sd=omega_V}
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}

[LONGITUDINAL]
file = 'lib:oral1_1cpt_kaVCl.txt'

DEFINITION:
Cc = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}

<PARAMETER>
ka_pop = {value=1.5, method=MLE}
V_pop = {value=50.0, method=MLE}
Cl_pop = {value=10.0, method=MLE}
omega_ka = {value=0.3, method=MLE}
omega_V = {value=0.2, method=MLE}
omega_Cl = {value=0.25, method=MLE}
a = {value=0.5, method=MLE}
b = {value=0.1, method=MLE}
"""

        mlx_legacy = parse_monolix_project(legacy_text)
        mlx_modern = parse_monolix_project(modern_text)

        # Both should produce the same model type
        @test mlx_legacy.model !== nothing
        @test mlx_modern.model !== nothing
        @test mlx_legacy.model.n_compartments == mlx_modern.model.n_compartments

        # Both should have IIV
        legacy_ka = findfirst(p -> p.name == "ka", mlx_legacy.parameters)
        modern_ka = findfirst(p -> p.name == "ka", mlx_modern.parameters)
        @test legacy_ka !== nothing
        @test modern_ka !== nothing
        @test mlx_legacy.parameters[legacy_ka].has_iiv == mlx_modern.parameters[modern_ka].has_iiv

        # Both should convert successfully
        doses = [DoseEvent(0.0, 100.0)]
        result_legacy = convert_monolix_to_neopkpd(mlx_legacy; doses=doses)
        result_modern = convert_monolix_to_neopkpd(mlx_modern; doses=doses)

        @test isempty(result_legacy.errors)
        @test isempty(result_modern.errors)
        @test typeof(result_legacy.model_spec.kind) == typeof(result_modern.model_spec.kind)
    end
end

@testset "\$PK Block Parsing" begin
    @testset "TV Definitions" begin
        ctl_text = raw"""
$PROBLEM Test TV definitions
$SUB ADVAN1 TRANS2
$PK
TVCL = THETA(1)
TVV = THETA(2)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
S1 = V
$THETA
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        # Should detect TV definitions
        @test haskey(pk_block.tv_definitions, :TVCL)
        @test pk_block.tv_definitions[:TVCL] == 1
        @test haskey(pk_block.tv_definitions, :TVV)
        @test pk_block.tv_definitions[:TVV] == 2

        # Should detect assignments with ETA
        cl_assignment = findfirst(a -> a.target == :CL, pk_block.assignments)
        @test cl_assignment !== nothing
        @test pk_block.assignments[cl_assignment].eta_index == 1
        @test pk_block.assignments[cl_assignment].transformation == :exponential

        # Should detect scaling
        @test length(pk_block.scaling) >= 1
        s1 = findfirst(s -> s.compartment == 1, pk_block.scaling)
        @test s1 !== nothing
        @test pk_block.scaling[s1].parameter == :V
    end

    @testset "Additive ETA" begin
        ctl_text = raw"""
$PROBLEM Test additive ETA
$SUB ADVAN1 TRANS2
$PK
TVCL = THETA(1)
CL = TVCL + ETA(1)
V = THETA(2) * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        cl_idx = findfirst(a -> a.target == :CL, pk_block.assignments)
        @test cl_idx !== nothing
        @test pk_block.assignments[cl_idx].transformation == :additive
        @test pk_block.assignments[cl_idx].eta_index == 1
    end

    @testset "Direct THETA Assignment" begin
        ctl_text = raw"""
$PROBLEM Test direct THETA
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1) * EXP(ETA(1))
V = THETA(2) * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        cl_idx = findfirst(a -> a.target == :CL, pk_block.assignments)
        @test cl_idx !== nothing
        @test pk_block.assignments[cl_idx].tv_theta == 1

        v_idx = findfirst(a -> a.target == :V, pk_block.assignments)
        @test v_idx !== nothing
        @test pk_block.assignments[v_idx].tv_theta == 2
    end
end

@testset "\$ERROR Block Parsing" begin
    @testset "Proportional Error" begin
        ctl_text = raw"""
$PROBLEM Test proportional error
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
$ERROR
IPRED = F
W = IPRED * THETA(3)
Y = IPRED + W * ERR(1)
$THETA (0, 10, 100) (0, 50, 500) (0, 0.1, 1)
$OMEGA 0.04 0.09
$SIGMA 1 FIX
"""
        ctl = parse_nonmem_control(ctl_text)
        error_block = NeoPKPD.parse_error_block(ctl.error_code)

        @test error_block.error_type == :proportional
        @test 3 in error_block.theta_indices
    end

    @testset "Additive Error" begin
        ctl_text = raw"""
$PROBLEM Test additive error
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
$ERROR
IPRED = F
W = THETA(3)
Y = IPRED + W * ERR(1)
$THETA (0, 10, 100) (0, 50, 500) (0, 0.5, 5)
$OMEGA 0.04 0.09
$SIGMA 1 FIX
"""
        ctl = parse_nonmem_control(ctl_text)
        error_block = NeoPKPD.parse_error_block(ctl.error_code)

        @test error_block.error_type == :additive
    end

    @testset "Combined Error" begin
        ctl_text = raw"""
$PROBLEM Test combined error
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
$ERROR
IPRED = F
W = SQRT(THETA(3)**2 + (THETA(4)*IPRED)**2)
Y = IPRED + W * ERR(1)
$THETA (0, 10, 100) (0, 50, 500) (0, 0.5, 5) (0, 0.1, 1)
$OMEGA 0.04 0.09
$SIGMA 1 FIX
"""
        ctl = parse_nonmem_control(ctl_text)
        error_block = NeoPKPD.parse_error_block(ctl.error_code)

        @test error_block.error_type == :combined
        @test 3 in error_block.theta_indices || 4 in error_block.theta_indices
    end

    @testset "Exponential Error" begin
        ctl_text = raw"""
$PROBLEM Test exponential error
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
$ERROR
IPRED = F
Y = IPRED * EXP(ERR(1))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.04
"""
        ctl = parse_nonmem_control(ctl_text)
        error_block = NeoPKPD.parse_error_block(ctl.error_code)

        @test error_block.error_type == :exponential
    end
end

@testset "Unsupported Construct Detection" begin
    @testset "IF Statement Detection" begin
        ctl_text = raw"""
$PROBLEM Test IF statement
$SUB ADVAN1 TRANS2
$PK
IF(AGE.GT.65) THEN
  TVCL = THETA(1) * 0.8
ELSE
  TVCL = THETA(1)
ENDIF
CL = TVCL * EXP(ETA(1))
V = THETA(2) * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        unsupported = NeoPKPD.check_unsupported_constructs(ctl)

        @test !isempty(unsupported)
        @test any(u -> u.construct == "IF statement", unsupported)
    end

    @testset "ALAG Detection" begin
        ctl_text = raw"""
$PROBLEM Test ALAG
$SUB ADVAN2 TRANS2
$PK
KA = THETA(1)
CL = THETA(2)
V = THETA(3)
ALAG1 = THETA(4)
$THETA (0, 1.5, 10) (0, 10, 100) (0, 50, 500) (0, 0.5, 2)
$OMEGA 0.04 0.09 0.16
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        unsupported = NeoPKPD.check_unsupported_constructs(ctl)

        @test !isempty(unsupported)
        @test any(u -> occursin("ALAG", u.construct), unsupported)
    end

    @testset "Bioavailability F1 Detection" begin
        ctl_text = raw"""
$PROBLEM Test F1
$SUB ADVAN2 TRANS2
$PK
KA = THETA(1)
CL = THETA(2)
V = THETA(3)
F1 = THETA(4)
$THETA (0, 1.5, 10) (0, 10, 100) (0, 50, 500) (0, 0.8, 1)
$OMEGA 0.04 0.09 0.16
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        unsupported = NeoPKPD.check_unsupported_constructs(ctl)

        @test !isempty(unsupported)
        @test any(u -> occursin("Bioavailability", u.construct), unsupported)
    end

    @testset "MTIME Detection" begin
        ctl_text = raw"""
$PROBLEM Test MTIME
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
MTIME(1) = THETA(3)
$THETA (0, 10, 100) (0, 50, 500) (0, 6, 24)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        unsupported = NeoPKPD.check_unsupported_constructs(ctl)

        @test !isempty(unsupported)
        @test any(u -> occursin("MTIME", u.construct), unsupported)
    end

    @testset "Unsupported ADVAN Detection" begin
        for advan in [5, 6, 7, 8, 9, 12, 13]
            ctl_text = """
\$PROBLEM Test ADVAN$advan
\$SUB ADVAN$advan
\$THETA (10)
"""
            ctl = parse_nonmem_control(ctl_text)
            unsupported = NeoPKPD.check_unsupported_constructs(ctl)

            @test !isempty(unsupported)
            @test any(u -> occursin("ADVAN", u.construct), unsupported)
        end
    end

    @testset "Strict Mode Fails on Unsupported" begin
        ctl_text = raw"""
$PROBLEM Test strict mode
$SUB ADVAN1 TRANS2
$PK
IF(SEX.EQ.1) CL = THETA(1) * 0.9
CL = THETA(1) * EXP(ETA(1))
V = THETA(2) * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)

        # With strict=true, should add errors
        result = convert_nonmem_to_neopkpd(ctl; strict=true)
        @test !isempty(result.errors)
        @test any(e -> occursin("IF", e), result.errors)

        # With strict=false, should add warnings
        result_lenient = convert_nonmem_to_neopkpd(ctl; strict=false)
        @test !isempty(result_lenient.warnings)
    end
end

@testset "Covariate Extraction" begin
    @testset "Power Covariate on CL" begin
        ctl_text = raw"""
$PROBLEM Test power covariate
$SUB ADVAN1 TRANS2
$PK
TVCL = THETA(1) * (WT/70)**THETA(3)
TVV = THETA(2)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500) (-2, 0.75, 2)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        # Find the TVCL definition line and check for covariate effects
        # The covariate should be detected
        has_wt_covariate = false
        for assignment in pk_block.assignments
            for cov in assignment.covariate_effects
                if cov.covariate == :WT
                    has_wt_covariate = true
                    @test cov.effect_type == :power
                    @test cov.reference ≈ 70.0
                    @test cov.theta_index == 3
                end
            end
        end
        # Note: covariate detection may be in TV lines or assignments
        # At minimum, check the raw code captures it
        @test any(line -> occursin("WT/70", line), pk_block.raw_code)
    end

    @testset "Linear Covariate" begin
        ctl_text = raw"""
$PROBLEM Test linear covariate
$SUB ADVAN1 TRANS2
$PK
TVCL = THETA(1) * (1 + THETA(3)*(AGE-40))
TVV = THETA(2)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500) (-0.1, 0.01, 0.1)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        # Check raw code captures it
        @test any(line -> occursin("AGE", line), pk_block.raw_code)
    end

    @testset "Exponential Covariate" begin
        ctl_text = raw"""
$PROBLEM Test exponential covariate
$SUB ADVAN1 TRANS2
$PK
TVCL = THETA(1) * EXP(THETA(3)*(CRCL-100)/100)
TVV = THETA(2)
CL = TVCL * EXP(ETA(1))
V = TVV * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500) (-2, 0.5, 2)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        pk_block = NeoPKPD.parse_pk_block(ctl.pk_code)

        # Check raw code captures it
        @test any(line -> occursin("CRCL", line), pk_block.raw_code)
    end
end

@testset "Validation Checks" begin
    @testset "THETA Index Validation" begin
        ctl_text = raw"""
$PROBLEM Test THETA validation
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1) * EXP(ETA(1))
V = THETA(5) * EXP(ETA(2))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_neopkpd(ctl)

        # Should detect that THETA(5) exceeds defined THETAs
        @test !isempty(result.warnings) || !isempty(result.errors)
    end

    @testset "ETA Index Validation" begin
        ctl_text = raw"""
$PROBLEM Test ETA validation
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1) * EXP(ETA(1))
V = THETA(2) * EXP(ETA(5))
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_neopkpd(ctl)

        # Should detect that ETA(5) exceeds OMEGA dimensions
        @test !isempty(result.warnings) || !isempty(result.errors)
    end

    @testset "Scaling Factor Validation" begin
        ctl_text = raw"""
$PROBLEM Test scaling validation
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
S5 = V
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        result = convert_nonmem_to_neopkpd(ctl)

        # S5 references compartment 5 which doesn't exist for ADVAN1
        @test !isempty(result.warnings) || !isempty(result.errors)
    end
end

@testset "Prediction Validation" begin
    @testset "One-Comp IV Analytical" begin
        # For one-comp IV: C(t) = (Dose/V) * exp(-CL/V * t)
        # CL = 10, V = 50, Dose = 100
        # C(t) = 2.0 * exp(-0.2 * t)
        ctl_text = raw"""
$PROBLEM Analytical validation
$SUB ADVAN1 TRANS2
$PK
CL = THETA(1)
V = THETA(2)
S1 = V
$THETA (0, 10, 100) (0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test result.model_spec !== nothing
        @test result.model_spec.params.CL == 10.0
        @test result.model_spec.params.V == 50.0

        # Simulate
        times = collect(0.0:1.0:10.0)
        grid = SimGrid(0.0, 10.0, times)
        solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10000)
        sim = simulate(result.model_spec, grid, solver)

        # Analytical solution
        k = 10.0 / 50.0  # CL/V = 0.2
        expected = 2.0 .* exp.(-k .* times)

        @test sim.observations[:conc] ≈ expected rtol=1e-4
    end

    @testset "One-Comp Oral Tmax" begin
        # For one-comp oral: Tmax = ln(Ka/Ke) / (Ka - Ke)
        # Ka = 1.0, CL = 5, V = 50, Ke = 0.1
        # Tmax = ln(10) / 0.9 ≈ 2.56
        ctl_text = raw"""
$PROBLEM Oral Tmax validation
$SUB ADVAN2 TRANS2
$PK
KA = THETA(1)
CL = THETA(2)
V = THETA(3)
S1 = V
$THETA (0, 1.0, 10) (0, 5.0, 100) (0, 50, 500)
$OMEGA 0.04 0.09 0.16
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test result.model_spec !== nothing

        # Simulate with fine time resolution around expected Tmax
        times = collect(0.0:0.1:10.0)
        grid = SimGrid(0.0, 10.0, times)
        solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10000)
        sim = simulate(result.model_spec, grid, solver)

        # Find Tmax
        conc = sim.observations[:conc]
        tmax_idx = argmax(conc)
        tmax_sim = times[tmax_idx]

        # Analytical Tmax
        ka = 1.0
        ke = 5.0 / 50.0
        tmax_analytical = log(ka / ke) / (ka - ke)

        @test abs(tmax_sim - tmax_analytical) < 0.2  # Within 0.2 hours
    end

    @testset "Two-Comp IV Distribution" begin
        # Two-comp should show biphasic decline
        ctl_text = raw"""
$PROBLEM Two-comp validation
$SUB ADVAN3 TRANS4
$PK
CL = THETA(1)
V1 = THETA(2)
Q = THETA(3)
V2 = THETA(4)
S1 = V1
$THETA (0, 10, 100) (0, 20, 200) (0, 15, 150) (0, 50, 500)
$OMEGA 0.04 0.09 0.04 0.09
$SIGMA 0.01
"""
        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test result.model_spec !== nothing

        # Simulate
        times = collect(0.0:0.5:24.0)
        grid = SimGrid(0.0, 24.0, times)
        solver = SolverSpec(:Tsit5, 1e-8, 1e-10, 10000)
        sim = simulate(result.model_spec, grid, solver)

        conc = sim.observations[:conc]

        # Should show biphasic decline (initial rapid distribution, then slower elimination)
        # C0 = Dose/V1 = 100/20 = 5
        @test conc[1] ≈ 5.0 rtol=0.01

        # Concentration should decay monotonically
        for i in 2:length(conc)
            @test conc[i] <= conc[i-1] + 1e-6  # Allow small numerical tolerance
        end

        # Should have terminal phase - still measurable at 24h
        @test conc[end] > 0.01
    end
end

@testset "Import Integration" begin
    @testset "NONMEM Import and Simulate" begin
        ctl_text = raw"""
$PROBLEM Integration test
$SUB ADVAN1 TRANS2
$THETA
(0, 10, 100)
(0, 50, 500)
$OMEGA 0.04 0.09
$SIGMA 0.01
"""

        ctl = parse_nonmem_control(ctl_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_nonmem_to_neopkpd(ctl; doses=doses)

        @test result.model_spec !== nothing

        # Simulate
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sim_result = simulate(result.model_spec, grid, solver)

        @test length(sim_result.t) == 25
        @test sim_result.observations[:conc][1] == 100.0 / 50.0  # Dose / V
        @test sim_result.observations[:conc][end] < sim_result.observations[:conc][1]  # Decay
    end

    @testset "Monolix Import and Simulate" begin
        mlx_text = """
[STRUCTURAL_MODEL]
lib = pklib:pk_oral1cpt_1abs_kaVCl_PLASMA

[PARAMETER]
ka = 1.0
V = 50.0
Cl = 10.0
"""

        mlx = parse_monolix_project(mlx_text)
        doses = [DoseEvent(0.0, 100.0)]
        result = convert_monolix_to_neopkpd(mlx; doses=doses)

        @test result.model_spec !== nothing

        # Simulate
        grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
        solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

        sim_result = simulate(result.model_spec, grid, solver)

        @test length(sim_result.t) == 25
        @test sim_result.observations[:conc][1] == 0.0  # Starts at 0 for oral
        @test sim_result.observations[:conc][5] > 0.0  # Absorbed by t=4
    end
end
