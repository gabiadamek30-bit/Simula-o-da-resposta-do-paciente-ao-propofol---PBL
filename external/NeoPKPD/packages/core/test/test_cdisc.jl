# Test suite for CDISC/SDTM Data Support

using Test
using NeoPKPD

@testset "CDISC Data Support" begin

    @testset "CDISC Types" begin
        # Test PCRecord
        pc1 = PCRecord(
            usubjid="SUBJ001",
            pctestcd="DRUG1",
            pcstresn=100.5,
            pctptnum=2.0
        )
        @test pc1.usubjid == "SUBJ001"
        @test pc1.pctestcd == "DRUG1"
        @test pc1.pcstresn == 100.5
        @test pc1.pctptnum == 2.0
        @test pc1.pcspec == "PLASMA"  # Default

        # Test EXRecord
        ex1 = EXRecord(
            usubjid="SUBJ001",
            exdose=100.0,
            exroute="INTRAVENOUS",
            exstdtc="2024-01-01T08:00:00"
        )
        @test ex1.usubjid == "SUBJ001"
        @test ex1.exdose == 100.0
        @test ex1.exroute == "INTRAVENOUS"
        @test ex1.exdosu == "mg"  # Default

        # Test DMRecord
        dm1 = DMRecord(
            usubjid="SUBJ001",
            age=45.0,
            sex="M",
            race="WHITE"
        )
        @test dm1.usubjid == "SUBJ001"
        @test dm1.age == 45.0
        @test dm1.sex == "M"
        @test dm1.race == "WHITE"
        @test dm1.ageu == "YEARS"  # Default

        # Test PPRecord
        pp1 = PPRecord(
            usubjid="SUBJ001",
            pptestcd="AUCLST",
            ppstresn=500.0,
            ppstresu="ng*h/mL"
        )
        @test pp1.usubjid == "SUBJ001"
        @test pp1.pptestcd == "AUCLST"
        @test pp1.ppstresn == 500.0

        # Test CDISCDataset
        dataset = CDISCDataset(
            pc=[pc1],
            ex=[ex1],
            dm=[dm1],
            pp=[pp1],
            study_id="STUDY001"
        )
        @test length(dataset.pc) == 1
        @test length(dataset.ex) == 1
        @test length(dataset.dm) == 1
        @test length(dataset.pp) == 1
        @test dataset.study_id == "STUDY001"

        # Test SubjectData
        doses = [DoseEvent(0.0, 100.0)]
        subj = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0, 4.0],
            [0.0, 50.0, 30.0, 15.0],
            doses;
            covariates=Dict{Symbol,Any}(:age => 45.0, :sex => "M"),
            lloq=0.5
        )
        @test subj.subject_id == "SUBJ001"
        @test length(subj.times) == 4
        @test length(subj.observations) == 4
        @test length(subj.doses) == 1
        @test subj.covariates[:age] == 45.0
        @test subj.lloq == 0.5

        # Test ObservedData
        obs = ObservedData(
            [subj];
            study_id="STUDY001",
            analyte="DRUG1",
            units="ng/mL",
            time_units="h"
        )
        @test length(obs.subjects) == 1
        @test obs.study_id == "STUDY001"
        @test obs.analyte == "DRUG1"
        @test obs.units == "ng/mL"
        @test obs.time_units == "h"
    end

    @testset "ObservedData Accessors" begin
        doses1 = [DoseEvent(0.0, 100.0)]
        doses2 = [DoseEvent(0.0, 200.0)]

        subj1 = SubjectData(
            "SUBJ001",
            [0.0, 1.0, 2.0],
            [0.0, 50.0, 30.0],
            doses1
        )
        subj2 = SubjectData(
            "SUBJ002",
            [0.0, 1.0, 2.0, 4.0],
            [0.0, 80.0, 50.0, 20.0],
            doses2
        )

        obs = ObservedData([subj1, subj2])

        @test subject_ids(obs) == ["SUBJ001", "SUBJ002"]
        @test n_subjects(obs) == 2
        @test n_observations(obs) == 7  # 3 + 4
        @test length(all_times(obs)) == 7
        @test length(all_observations(obs)) == 7
    end

    @testset "ISO 8601 Datetime Parsing" begin
        # Full datetime
        t1 = parse_iso8601_datetime("2024-01-15T10:30:45")
        @test t1 !== nothing
        @test t1 > 0

        # Date only
        t2 = parse_iso8601_datetime("2024-01-15")
        @test t2 !== nothing
        @test t2 > 0

        # Empty/invalid
        @test parse_iso8601_datetime("") === nothing
        @test parse_iso8601_datetime(".") === nothing
        @test parse_iso8601_datetime("NA") === nothing
        @test parse_iso8601_datetime("invalid") === nothing
    end

    @testset "ISO 8601 Duration Parsing" begin
        # Hours
        @test parse_iso8601_duration("PT1H") == 1.0
        @test parse_iso8601_duration("PT2H") == 2.0
        @test parse_iso8601_duration("PT1.5H") == 1.5

        # Minutes
        @test parse_iso8601_duration("PT30M") == 0.5
        @test parse_iso8601_duration("PT60M") == 1.0

        # Combined
        @test parse_iso8601_duration("PT1H30M") == 1.5

        # Days
        @test parse_iso8601_duration("P1D") == 24.0
        @test parse_iso8601_duration("P2D") == 48.0

        # Seconds
        @test parse_iso8601_duration("PT3600S") == 1.0

        # Empty/invalid
        @test parse_iso8601_duration("") == 0.0
        @test parse_iso8601_duration(".") == 0.0
        @test parse_iso8601_duration("NA") == 0.0
    end

    @testset "CDISCConversionConfig" begin
        # Default config
        config1 = CDISCConversionConfig()
        @test config1.time_reference == :first_dose
        @test config1.time_units == :hours
        @test config1.include_blq == true
        @test config1.blq_handling == :half_lloq
        @test config1.analyte_filter == ""

        # Custom config
        config2 = CDISCConversionConfig(
            time_reference=:reference_time,
            time_units=:days,
            include_blq=false,
            blq_handling=:exclude
        )
        @test config2.time_reference == :reference_time
        @test config2.time_units == :days
        @test config2.include_blq == false
        @test config2.blq_handling == :exclude
    end

    @testset "CDISC to NeoPKPD Conversion" begin
        # Create test CDISC dataset
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0, pcstresu="ng/mL"),
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=50.0, pctptnum=2.0, pcstresu="ng/mL"),
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=25.0, pctptnum=4.0, pcstresu="ng/mL"),
            PCRecord(usubjid="SUBJ002", pctestcd="DRUG1", pcstresn=120.0, pctptnum=0.0, pcstresu="ng/mL"),
            PCRecord(usubjid="SUBJ002", pctestcd="DRUG1", pcstresn=60.0, pctptnum=2.0, pcstresu="ng/mL"),
        ]

        ex_records = [
            EXRecord(usubjid="SUBJ001", exdose=100.0, exdy=1),
            EXRecord(usubjid="SUBJ002", exdose=200.0, exdy=1),
        ]

        dm_records = [
            DMRecord(usubjid="SUBJ001", age=35.0, sex="M", race="WHITE"),
            DMRecord(usubjid="SUBJ002", age=45.0, sex="F", race="ASIAN"),
        ]

        dataset = CDISCDataset(
            pc=pc_records,
            ex=ex_records,
            dm=dm_records,
            study_id="TEST001"
        )

        # Convert to NeoPKPD format
        result = cdisc_to_observed(dataset)

        @test isempty(result.errors)
        @test result.observed_data !== nothing
        @test result.n_subjects == 2
        @test result.n_observations == 5

        obs = result.observed_data
        @test obs.study_id == "TEST001"
        @test obs.analyte == "DRUG1"
        @test obs.units == "ng/mL"
    end

    @testset "CDISC Conversion - Empty Dataset" begin
        dataset = CDISCDataset()
        result = cdisc_to_observed(dataset)

        @test !isempty(result.errors)
        @test result.observed_data === nothing
        @test result.n_subjects == 0
    end

    @testset "CDISC Conversion - Analyte Filter" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0),
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG2", pcstresn=50.0, pctptnum=0.0),
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=25.0, pctptnum=2.0),
        ]

        dataset = CDISCDataset(pc=pc_records, study_id="TEST002")

        # Filter to DRUG1 only
        config = CDISCConversionConfig(analyte_filter="DRUG1")
        result = cdisc_to_observed(dataset; config=config)

        @test result.n_observations == 2
        @test result.observed_data.analyte == "DRUG1"
    end

    @testset "CDISC Conversion - BLQ Handling" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0, pclloq=1.0),
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=nothing, pctptnum=12.0, pcstat="BLQ", pclloq=1.0),
        ]

        dataset = CDISCDataset(pc=pc_records)

        # Half LLOQ (default)
        config1 = CDISCConversionConfig(blq_handling=:half_lloq)
        result1 = cdisc_to_observed(dataset; config=config1)
        @test result1.n_observations == 2
        @test result1.observed_data.subjects[1].observations[2] == 0.5  # half of LLOQ=1.0

        # Zero
        config2 = CDISCConversionConfig(blq_handling=:zero)
        result2 = cdisc_to_observed(dataset; config=config2)
        @test result2.observed_data.subjects[1].observations[2] == 0.0

        # Exclude
        config3 = CDISCConversionConfig(include_blq=false, blq_handling=:exclude)
        result3 = cdisc_to_observed(dataset; config=config3)
        @test result3.n_observations == 1
        @test result3.n_excluded == 1
    end

    @testset "CDISC Conversion - Time Units" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=24.0),  # 24 hours
        ]

        dataset = CDISCDataset(pc=pc_records)

        # Hours (default)
        config1 = CDISCConversionConfig(time_units=:hours)
        result1 = cdisc_to_observed(dataset; config=config1)
        @test result1.observed_data.subjects[1].times[1] == 24.0
        @test result1.observed_data.time_units == "h"

        # Days
        config2 = CDISCConversionConfig(time_units=:days)
        result2 = cdisc_to_observed(dataset; config=config2)
        @test result2.observed_data.subjects[1].times[1] == 1.0  # 24h = 1 day
        @test result2.observed_data.time_units == "d"

        # Minutes
        config3 = CDISCConversionConfig(time_units=:minutes)
        result3 = cdisc_to_observed(dataset; config=config3)
        @test result3.observed_data.subjects[1].times[1] == 1440.0  # 24h = 1440 min
        @test result3.observed_data.time_units == "min"
    end

    @testset "CDISC Conversion - Covariate Extraction" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0),
        ]

        dm_records = [
            DMRecord(
                usubjid="SUBJ001",
                age=45.0,
                sex="M",
                race="WHITE",
                ethnic="NOT HISPANIC",
                armcd="TRT01",
                country="USA"
            ),
        ]

        dataset = CDISCDataset(pc=pc_records, dm=dm_records)
        result = cdisc_to_observed(dataset)

        @test result.observed_data !== nothing
        covs = result.observed_data.subjects[1].covariates

        @test covs[:age] == 45.0
        @test covs[:sex] == "M"
        @test covs[:race] == "WHITE"
        @test covs[:ethnic] == "NOT HISPANIC"
        @test covs[:arm] == "TRT01"
        @test covs[:country] == "USA"
    end

    @testset "CDISC Conversion - Dose Extraction" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0),
        ]

        ex_records = [
            EXRecord(usubjid="SUBJ001", exdose=100.0, exstdtc="2024-01-01T08:00:00"),
            EXRecord(usubjid="SUBJ001", exdose=100.0, exstdtc="2024-01-02T08:00:00"),
        ]

        dataset = CDISCDataset(pc=pc_records, ex=ex_records)
        result = cdisc_to_observed(dataset)

        @test result.observed_data !== nothing
        doses = result.observed_data.subjects[1].doses
        @test length(doses) == 2
        @test doses[1].amount == 100.0
        @test doses[2].amount == 100.0
    end

    @testset "CDISC Conversion - Infusion Duration" begin
        pc_records = [
            PCRecord(usubjid="SUBJ001", pctestcd="DRUG1", pcstresn=100.0, pctptnum=0.0),
        ]

        ex_records = [
            EXRecord(usubjid="SUBJ001", exdose=100.0, exstdtc="2024-01-01T08:00:00", exdur="PT1H"),
        ]

        dataset = CDISCDataset(pc=pc_records, ex=ex_records)
        result = cdisc_to_observed(dataset)

        @test result.observed_data !== nothing
        doses = result.observed_data.subjects[1].doses
        @test length(doses) == 1
        @test doses[1].amount == 100.0
        @test doses[1].duration == 1.0  # 1 hour infusion
    end

end
