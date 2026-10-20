# CDISC/SDTM Data Types
# Structures for representing CDISC standard domains (PC, EX, DM, PP)

export PCRecord, EXRecord, DMRecord, PPRecord
export CDISCDataset, ObservedData, SubjectData

"""
PC Domain Record - Pharmacokinetic Concentrations.

Standard CDISC SDTM PC domain variables for drug concentration measurements.

Fields:
- studyid: Study identifier
- domain: Domain abbreviation (always "PC")
- usubjid: Unique subject identifier
- pcseq: Sequence number
- pcgrpid: Group ID (optional)
- pcrefid: Reference ID (optional)
- pcspid: Sponsor-defined ID (optional)
- pctestcd: Test short name (e.g., "DRUG1")
- pctest: Test name (e.g., "Drug 1 Concentration")
- pccat: Category (optional)
- pcorres: Result or finding as originally collected
- pcorresu: Original units
- pcstresc: Character result/finding in standard format
- pcstresn: Numeric result/finding in standard units
- pcstresu: Standard units
- pcstat: Completion status
- pcreasnd: Reason not done
- pcnam: Vendor name (optional)
- pcspec: Specimen material type
- pcspccnd: Specimen condition (optional)
- pcmethod: Method of test (optional)
- pcblfl: Baseline flag
- pcdrvfl: Derived flag (optional)
- pclloq: Lower limit of quantitation
- pcdy: Study day
- pctpt: Planned time point name
- pctptnum: Planned time point number
- pceltm: Planned elapsed time from reference
- pctptref: Time point reference
- pcrftdtc: Date/time of reference point (ISO 8601)
- pcstdtc: Start date/time of observation (ISO 8601)
"""
struct PCRecord
    studyid::String
    usubjid::String
    pcseq::Int
    pctestcd::String
    pctest::String
    pcorres::String
    pcorresu::String
    pcstresn::Union{Nothing,Float64}
    pcstresu::String
    pcstat::String
    pclloq::Union{Nothing,Float64}
    pcspec::String
    pcdy::Union{Nothing,Int}
    pctpt::String
    pctptnum::Union{Nothing,Float64}
    pceltm::String
    pcrftdtc::String
    pcstdtc::String
    pcblfl::String
    pcmethod::String

    function PCRecord(;
        studyid::String="",
        usubjid::String="",
        pcseq::Int=0,
        pctestcd::String="",
        pctest::String="",
        pcorres::String="",
        pcorresu::String="",
        pcstresn::Union{Nothing,Float64}=nothing,
        pcstresu::String="",
        pcstat::String="",
        pclloq::Union{Nothing,Float64}=nothing,
        pcspec::String="PLASMA",
        pcdy::Union{Nothing,Int}=nothing,
        pctpt::String="",
        pctptnum::Union{Nothing,Float64}=nothing,
        pceltm::String="",
        pcrftdtc::String="",
        pcstdtc::String="",
        pcblfl::String="",
        pcmethod::String=""
    )
        new(studyid, usubjid, pcseq, pctestcd, pctest, pcorres, pcorresu,
            pcstresn, pcstresu, pcstat, pclloq, pcspec, pcdy, pctpt, pctptnum,
            pceltm, pcrftdtc, pcstdtc, pcblfl, pcmethod)
    end
end

"""
EX Domain Record - Exposure.

Standard CDISC SDTM EX domain variables for drug administration/dosing.

Fields:
- studyid: Study identifier
- usubjid: Unique subject identifier
- exseq: Sequence number
- extrt: Name of treatment
- excat: Category of treatment (optional)
- exdose: Dose
- exdosu: Dose units
- exdosfrm: Dose form
- exdosfrq: Dosing frequency per interval
- exroute: Route of administration
- exstdtc: Start date/time (ISO 8601)
- exendtc: End date/time (ISO 8601)
- exdy: Study day of start
- exendy: Study day of end
- exdur: Duration (ISO 8601 duration format)
"""
struct EXRecord
    studyid::String
    usubjid::String
    exseq::Int
    extrt::String
    excat::String
    exdose::Float64
    exdosu::String
    exdosfrm::String
    exdosfrq::String
    exroute::String
    exstdtc::String
    exendtc::String
    exdy::Union{Nothing,Int}
    exendy::Union{Nothing,Int}
    exdur::String

    function EXRecord(;
        studyid::String="",
        usubjid::String="",
        exseq::Int=0,
        extrt::String="",
        excat::String="",
        exdose::Float64=0.0,
        exdosu::String="mg",
        exdosfrm::String="",
        exdosfrq::String="",
        exroute::String="",
        exstdtc::String="",
        exendtc::String="",
        exdy::Union{Nothing,Int}=nothing,
        exendy::Union{Nothing,Int}=nothing,
        exdur::String=""
    )
        new(studyid, usubjid, exseq, extrt, excat, exdose, exdosu, exdosfrm,
            exdosfrq, exroute, exstdtc, exendtc, exdy, exendy, exdur)
    end
end

"""
DM Domain Record - Demographics.

Standard CDISC SDTM DM domain variables for subject demographics/covariates.

Fields:
- studyid: Study identifier
- usubjid: Unique subject identifier
- subjid: Subject identifier for the study
- rfstdtc: Subject reference start date/time (ISO 8601)
- rfendtc: Subject reference end date/time (ISO 8601)
- siteid: Study site identifier
- brthdtc: Date/time of birth (ISO 8601)
- age: Age
- ageu: Age units
- sex: Sex (M/F/U)
- race: Race
- ethnic: Ethnicity
- armcd: Arm code
- arm: Description of arm
- country: Country
- dmdtc: Date/time of collection (ISO 8601)
- dmdy: Study day of collection
"""
struct DMRecord
    studyid::String
    usubjid::String
    subjid::String
    rfstdtc::String
    rfendtc::String
    siteid::String
    brthdtc::String
    age::Union{Nothing,Float64}
    ageu::String
    sex::String
    race::String
    ethnic::String
    armcd::String
    arm::String
    country::String
    dmdtc::String
    dmdy::Union{Nothing,Int}

    function DMRecord(;
        studyid::String="",
        usubjid::String="",
        subjid::String="",
        rfstdtc::String="",
        rfendtc::String="",
        siteid::String="",
        brthdtc::String="",
        age::Union{Nothing,Float64}=nothing,
        ageu::String="YEARS",
        sex::String="",
        race::String="",
        ethnic::String="",
        armcd::String="",
        arm::String="",
        country::String="",
        dmdtc::String="",
        dmdy::Union{Nothing,Int}=nothing
    )
        new(studyid, usubjid, subjid, rfstdtc, rfendtc, siteid, brthdtc,
            age, ageu, sex, race, ethnic, armcd, arm, country, dmdtc, dmdy)
    end
end

"""
PP Domain Record - Pharmacokinetic Parameters.

Standard CDISC SDTM PP domain for derived PK parameters (AUC, Cmax, etc.).

Fields:
- studyid: Study identifier
- usubjid: Unique subject identifier
- ppseq: Sequence number
- ppgrpid: Group ID
- pptestcd: PK parameter short name (e.g., "AUCLST", "CMAX")
- pptest: PK parameter name
- ppcat: Category for PK parameter
- pporres: Result as originally collected
- pporresu: Original units
- ppstresc: Character result in standard format
- ppstresn: Numeric result in standard units
- ppstresu: Standard units
- ppstat: Completion status
- ppreasnd: Reason not done
- ppspec: Specimen material type
- ppdtc: Date/time of parameter derivation
"""
struct PPRecord
    studyid::String
    usubjid::String
    ppseq::Int
    ppgrpid::String
    pptestcd::String
    pptest::String
    ppcat::String
    pporres::String
    pporresu::String
    ppstresn::Union{Nothing,Float64}
    ppstresu::String
    ppstat::String
    ppreasnd::String
    ppspec::String
    ppdtc::String

    function PPRecord(;
        studyid::String="",
        usubjid::String="",
        ppseq::Int=0,
        ppgrpid::String="",
        pptestcd::String="",
        pptest::String="",
        ppcat::String="",
        pporres::String="",
        pporresu::String="",
        ppstresn::Union{Nothing,Float64}=nothing,
        ppstresu::String="",
        ppstat::String="",
        ppreasnd::String="",
        ppspec::String="PLASMA",
        ppdtc::String=""
    )
        new(studyid, usubjid, ppseq, ppgrpid, pptestcd, pptest, ppcat,
            pporres, pporresu, ppstresn, ppstresu, ppstat, ppreasnd, ppspec, ppdtc)
    end
end

"""
Complete CDISC Dataset combining multiple domains.

Fields:
- pc: Pharmacokinetic Concentrations domain
- ex: Exposure domain
- dm: Demographics domain
- pp: Pharmacokinetic Parameters domain (optional, derived)
- study_id: Study identifier
"""
struct CDISCDataset
    pc::Vector{PCRecord}
    ex::Vector{EXRecord}
    dm::Vector{DMRecord}
    pp::Vector{PPRecord}
    study_id::String

    function CDISCDataset(;
        pc::Vector{PCRecord}=PCRecord[],
        ex::Vector{EXRecord}=EXRecord[],
        dm::Vector{DMRecord}=DMRecord[],
        pp::Vector{PPRecord}=PPRecord[],
        study_id::String=""
    )
        new(pc, ex, dm, pp, study_id)
    end
end

"""
Individual subject data extracted from CDISC domains.

Fields:
- subject_id: Unique subject identifier
- times: Time points (hours from first dose)
- observations: Observed concentrations
- doses: Vector of DoseEvent
- covariates: Subject-level covariates (age, weight, sex, etc.)
- lloq: Lower limit of quantitation
- blq_flags: Below LOQ flags for each observation
"""
struct SubjectData
    subject_id::String
    times::Vector{Float64}
    observations::Vector{Float64}
    doses::Vector{DoseEvent}
    covariates::Dict{Symbol,Any}
    lloq::Float64
    blq_flags::Vector{Bool}

    function SubjectData(
        subject_id::String,
        times::Vector{Float64},
        observations::Vector{Float64},
        doses::Vector{DoseEvent};
        covariates::Dict{Symbol,Any}=Dict{Symbol,Any}(),
        lloq::Float64=0.0,
        blq_flags::Vector{Bool}=Bool[]
    )
        if isempty(blq_flags)
            blq_flags = fill(false, length(observations))
        end
        new(subject_id, times, observations, doses, covariates, lloq, blq_flags)
    end
end

"""
Observed data structure for population analysis.

Contains concentration-time data from multiple subjects in a format
suitable for population modeling and VPC analysis.

Fields:
- subjects: Vector of SubjectData for each subject
- study_id: Study identifier
- analyte: Analyte/drug name
- units: Concentration units
- time_units: Time units
"""
struct ObservedData
    subjects::Vector{SubjectData}
    study_id::String
    analyte::String
    units::String
    time_units::String

    function ObservedData(
        subjects::Vector{SubjectData};
        study_id::String="",
        analyte::String="",
        units::String="ng/mL",
        time_units::String="h"
    )
        new(subjects, study_id, analyte, units, time_units)
    end
end

# Convenience accessors for ObservedData
"""Get all unique subject IDs from observed data."""
function subject_ids(data::ObservedData)::Vector{String}
    return [s.subject_id for s in data.subjects]
end

"""Get total number of subjects."""
function n_subjects(data::ObservedData)::Int
    return length(data.subjects)
end

"""Get total number of observations across all subjects."""
function n_observations(data::ObservedData)::Int
    return sum(length(s.observations) for s in data.subjects)
end

"""Get all times pooled across subjects."""
function all_times(data::ObservedData)::Vector{Float64}
    times = Float64[]
    for s in data.subjects
        append!(times, s.times)
    end
    return times
end

"""Get all observations pooled across subjects."""
function all_observations(data::ObservedData)::Vector{Float64}
    obs = Float64[]
    for s in data.subjects
        append!(obs, s.observations)
    end
    return obs
end

export subject_ids, n_subjects, n_observations, all_times, all_observations
