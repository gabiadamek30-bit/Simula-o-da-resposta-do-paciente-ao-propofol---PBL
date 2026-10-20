# ADaM (Analysis Data Model) Types
# CDISC ADaM standard datasets for analysis-ready data

export ADSLRecord, ADPCRecord, ADPPRecord, ADEXRecord
export ADaMDataset, ADaMPKDataset

# ============================================================================
# ADSL - Subject Level Analysis Dataset
# ============================================================================

"""
    ADSLRecord

ADSL (Subject Level Analysis Dataset) record.

Contains one record per subject with baseline and analysis covariates.
This is the foundational ADaM dataset that all other ADaM datasets can merge to.

Key Variables (CDISC ADaM IG):
- STUDYID: Study identifier
- USUBJID: Unique subject identifier
- SUBJID: Subject identifier for the study
- SITEID: Study site identifier
- AGE: Age (numeric)
- AGEU: Age units
- AGEGR1: Age group 1 (categorical)
- SEX: Sex (M/F/U)
- SEXN: Sex (numeric: 1=M, 2=F)
- RACE: Race
- RACEN: Race (numeric)
- ETHNIC: Ethnicity
- SAFFL: Safety population flag (Y/N)
- PKFL: PK population flag (Y/N)
- ITTFL: Intent-to-treat population flag (Y/N)
- TRT01P: Planned treatment for period 1
- TRT01A: Actual treatment for period 1
- TRTSDT: Treatment start date (numeric)
- TRTEDT: Treatment end date (numeric)
- Baseline covariates: BWEIGHT, BHEIGHT, BBMI, BCRCL, etc.
"""
struct ADSLRecord
    # Identifiers
    studyid::String
    usubjid::String
    subjid::String
    siteid::String

    # Demographics
    age::Union{Nothing, Float64}
    ageu::String
    agegr1::String
    sex::String
    sexn::Union{Nothing, Int}
    race::String
    racen::Union{Nothing, Int}
    ethnic::String

    # Population flags
    saffl::String     # Safety population
    pkfl::String      # PK population
    ittfl::String     # Intent-to-treat
    pkpopfl::String   # PK evaluable population

    # Treatment
    trt01p::String    # Planned treatment period 1
    trt01a::String    # Actual treatment period 1
    trt01pn::Union{Nothing, Int}   # Planned treatment (numeric)
    trt01an::Union{Nothing, Int}   # Actual treatment (numeric)
    trtsdt::Union{Nothing, Int}    # Treatment start date (SAS date)
    trtedt::Union{Nothing, Int}    # Treatment end date (SAS date)

    # Baseline covariates
    bweight::Union{Nothing, Float64}   # Baseline weight (kg)
    bheight::Union{Nothing, Float64}   # Baseline height (cm)
    bbmi::Union{Nothing, Float64}      # Baseline BMI
    bcrcl::Union{Nothing, Float64}     # Baseline creatinine clearance
    begfr::Union{Nothing, Float64}     # Baseline eGFR
    bast::Union{Nothing, Float64}      # Baseline AST
    balt::Union{Nothing, Float64}      # Baseline ALT
    bbili::Union{Nothing, Float64}     # Baseline bilirubin
    balb::Union{Nothing, Float64}      # Baseline albumin

    # Additional flags
    complfl::String   # Completer flag
    dcsreas::String   # Discontinuation reason

    function ADSLRecord(;
        studyid::String="",
        usubjid::String="",
        subjid::String="",
        siteid::String="",
        age::Union{Nothing, Float64}=nothing,
        ageu::String="YEARS",
        agegr1::String="",
        sex::String="",
        sexn::Union{Nothing, Int}=nothing,
        race::String="",
        racen::Union{Nothing, Int}=nothing,
        ethnic::String="",
        saffl::String="Y",
        pkfl::String="Y",
        ittfl::String="Y",
        pkpopfl::String="Y",
        trt01p::String="",
        trt01a::String="",
        trt01pn::Union{Nothing, Int}=nothing,
        trt01an::Union{Nothing, Int}=nothing,
        trtsdt::Union{Nothing, Int}=nothing,
        trtedt::Union{Nothing, Int}=nothing,
        bweight::Union{Nothing, Float64}=nothing,
        bheight::Union{Nothing, Float64}=nothing,
        bbmi::Union{Nothing, Float64}=nothing,
        bcrcl::Union{Nothing, Float64}=nothing,
        begfr::Union{Nothing, Float64}=nothing,
        bast::Union{Nothing, Float64}=nothing,
        balt::Union{Nothing, Float64}=nothing,
        bbili::Union{Nothing, Float64}=nothing,
        balb::Union{Nothing, Float64}=nothing,
        complfl::String="",
        dcsreas::String=""
    )
        new(studyid, usubjid, subjid, siteid, age, ageu, agegr1, sex, sexn,
            race, racen, ethnic, saffl, pkfl, ittfl, pkpopfl,
            trt01p, trt01a, trt01pn, trt01an, trtsdt, trtedt,
            bweight, bheight, bbmi, bcrcl, begfr, bast, balt, bbili, balb,
            complfl, dcsreas)
    end
end

# ============================================================================
# ADPC - PK Concentration Analysis Dataset
# ============================================================================

"""
    ADPCRecord

ADPC (PK Concentration Analysis Dataset) record.

Analysis-ready concentration data with derived variables for PK analysis.
Each record represents one concentration measurement.

Key Variables:
- STUDYID, USUBJID: Identifiers
- PARAM/PARAMCD: Parameter description/code (e.g., "Drug X Plasma Concentration")
- AVAL: Analysis value (concentration)
- AVALU: Analysis value units
- AVALC: Character analysis value
- ADT/ATM: Analysis date/time
- ATPT/ATPTN: Analysis timepoint (name/number)
- ARELTM: Actual relative time from reference
- NRELTM: Nominal relative time from reference
- AFRLT: Actual relative time from first dose
- NFRLT: Nominal relative time from first dose
- ARRLT: Actual relative time from most recent dose
- NRRLT: Nominal relative time from most recent dose
- DOSESSION: Dose session number
- PCTPT: PC timepoint name (from SDTM)
- ANL01FL: Analysis record flag 01
- BLQ: Below limit of quantitation flag
- LLOQ: Lower limit of quantitation
"""
struct ADPCRecord
    # Identifiers
    studyid::String
    usubjid::String
    pcseq::Int

    # Parameter identification
    paramcd::String   # Parameter code (e.g., "DRUGX")
    param::String     # Parameter description
    paramn::Union{Nothing, Int}

    # Analysis values
    aval::Union{Nothing, Float64}     # Analysis value (concentration)
    avalu::String                      # Analysis value units
    avalc::String                      # Character analysis value (for BLQ)
    base::Union{Nothing, Float64}     # Baseline value
    chg::Union{Nothing, Float64}      # Change from baseline
    pchg::Union{Nothing, Float64}     # Percent change from baseline

    # Time variables
    adt::Union{Nothing, Int}          # Analysis date (SAS date)
    atm::Union{Nothing, Float64}      # Analysis time (SAS time)
    adtm::Union{Nothing, Float64}     # Analysis datetime
    atpt::String                      # Analysis timepoint name
    atptn::Union{Nothing, Float64}    # Analysis timepoint number

    # Relative times (key PK variables)
    areltm::Union{Nothing, Float64}   # Actual relative time from reference
    nreltm::Union{Nothing, Float64}   # Nominal relative time from reference
    afrlt::Union{Nothing, Float64}    # Actual time from first dose
    nfrlt::Union{Nothing, Float64}    # Nominal time from first dose
    arrlt::Union{Nothing, Float64}    # Actual time from most recent dose
    nrrlt::Union{Nothing, Float64}    # Nominal time from most recent dose

    # Dosing information
    dosession::Union{Nothing, Int}    # Dose session number (occasion)
    doseno::Union{Nothing, Int}       # Dose number
    dose::Union{Nothing, Float64}     # Dose amount
    doseu::String                      # Dose units

    # Treatment
    trtp::String      # Planned treatment
    trta::String      # Actual treatment

    # Flags
    anl01fl::String   # Analysis record flag 01 (Y = include in primary analysis)
    anl02fl::String   # Analysis record flag 02
    blqfl::String     # Below LOQ flag (Y/N)
    ablfl::String     # Baseline record flag
    avisit::String    # Analysis visit
    avisitn::Union{Nothing, Int}

    # BLQ handling
    lloq::Union{Nothing, Float64}     # Lower limit of quantitation
    blq_handling::String              # How BLQ was handled (e.g., "LLOQ/2", "0", "MISSING")

    # Specimen
    spec::String      # Specimen type

    function ADPCRecord(;
        studyid::String="",
        usubjid::String="",
        pcseq::Int=0,
        paramcd::String="",
        param::String="",
        paramn::Union{Nothing, Int}=nothing,
        aval::Union{Nothing, Float64}=nothing,
        avalu::String="ng/mL",
        avalc::String="",
        base::Union{Nothing, Float64}=nothing,
        chg::Union{Nothing, Float64}=nothing,
        pchg::Union{Nothing, Float64}=nothing,
        adt::Union{Nothing, Int}=nothing,
        atm::Union{Nothing, Float64}=nothing,
        adtm::Union{Nothing, Float64}=nothing,
        atpt::String="",
        atptn::Union{Nothing, Float64}=nothing,
        areltm::Union{Nothing, Float64}=nothing,
        nreltm::Union{Nothing, Float64}=nothing,
        afrlt::Union{Nothing, Float64}=nothing,
        nfrlt::Union{Nothing, Float64}=nothing,
        arrlt::Union{Nothing, Float64}=nothing,
        nrrlt::Union{Nothing, Float64}=nothing,
        dosession::Union{Nothing, Int}=nothing,
        doseno::Union{Nothing, Int}=nothing,
        dose::Union{Nothing, Float64}=nothing,
        doseu::String="mg",
        trtp::String="",
        trta::String="",
        anl01fl::String="Y",
        anl02fl::String="",
        blqfl::String="N",
        ablfl::String="",
        avisit::String="",
        avisitn::Union{Nothing, Int}=nothing,
        lloq::Union{Nothing, Float64}=nothing,
        blq_handling::String="",
        spec::String="PLASMA"
    )
        new(studyid, usubjid, pcseq, paramcd, param, paramn,
            aval, avalu, avalc, base, chg, pchg,
            adt, atm, adtm, atpt, atptn,
            areltm, nreltm, afrlt, nfrlt, arrlt, nrrlt,
            dosession, doseno, dose, doseu,
            trtp, trta, anl01fl, anl02fl, blqfl, ablfl, avisit, avisitn,
            lloq, blq_handling, spec)
    end
end

# ============================================================================
# ADPP - PK Parameters Analysis Dataset
# ============================================================================

"""
    ADPPRecord

ADPP (PK Parameters Analysis Dataset) record.

Analysis-ready NCA/compartmental PK parameters.
Each record represents one PK parameter for one subject/treatment.

Key Variables:
- PARAM/PARAMCD: Parameter code (e.g., "AUCLST", "CMAX", "TMAX", "THALF")
- AVAL: Analysis value (parameter value)
- AVALU: Units
- ANL01FL: Analysis record flag
"""
struct ADPPRecord
    # Identifiers
    studyid::String
    usubjid::String
    ppseq::Int

    # Parameter
    paramcd::String   # Parameter code (AUCLST, CMAX, TMAX, etc.)
    param::String     # Full parameter name
    paramn::Union{Nothing, Int}
    parcat1::String   # Parameter category 1 (e.g., "NCA", "COMPARTMENTAL")
    parcat2::String   # Parameter category 2 (e.g., "PRIMARY", "SECONDARY")

    # Analysis values
    aval::Union{Nothing, Float64}     # Analysis value
    avalu::String                      # Units
    avalc::String                      # Character value
    base::Union{Nothing, Float64}     # Baseline (for bioequivalence)
    chg::Union{Nothing, Float64}      # Change from baseline

    # Treatment
    trtp::String
    trta::String
    trtpn::Union{Nothing, Int}
    trtan::Union{Nothing, Int}

    # Period/Sequence (for crossover studies)
    aperiod::Union{Nothing, Int}      # Analysis period
    aperiodc::String                   # Analysis period (character)
    aseq::Union{Nothing, Int}         # Analysis sequence

    # Visit
    avisit::String
    avisitn::Union{Nothing, Int}

    # Flags
    anl01fl::String   # Primary analysis flag
    anl02fl::String   # Secondary analysis flag
    ablfl::String     # Baseline flag

    # Dosing context
    doseno::Union{Nothing, Int}
    dose::Union{Nothing, Float64}
    doseu::String

    # Method
    method::String    # Method used (e.g., "LINEAR", "LOG-LINEAR")

    function ADPPRecord(;
        studyid::String="",
        usubjid::String="",
        ppseq::Int=0,
        paramcd::String="",
        param::String="",
        paramn::Union{Nothing, Int}=nothing,
        parcat1::String="NCA",
        parcat2::String="",
        aval::Union{Nothing, Float64}=nothing,
        avalu::String="",
        avalc::String="",
        base::Union{Nothing, Float64}=nothing,
        chg::Union{Nothing, Float64}=nothing,
        trtp::String="",
        trta::String="",
        trtpn::Union{Nothing, Int}=nothing,
        trtan::Union{Nothing, Int}=nothing,
        aperiod::Union{Nothing, Int}=nothing,
        aperiodc::String="",
        aseq::Union{Nothing, Int}=nothing,
        avisit::String="",
        avisitn::Union{Nothing, Int}=nothing,
        anl01fl::String="Y",
        anl02fl::String="",
        ablfl::String="",
        doseno::Union{Nothing, Int}=nothing,
        dose::Union{Nothing, Float64}=nothing,
        doseu::String="mg",
        method::String=""
    )
        new(studyid, usubjid, ppseq, paramcd, param, paramn, parcat1, parcat2,
            aval, avalu, avalc, base, chg,
            trtp, trta, trtpn, trtan,
            aperiod, aperiodc, aseq, avisit, avisitn,
            anl01fl, anl02fl, ablfl,
            doseno, dose, doseu, method)
    end
end

# ============================================================================
# ADEX - Exposure Analysis Dataset
# ============================================================================

"""
    ADEXRecord

ADEX (Exposure Analysis Dataset) record.

Analysis-ready exposure/dosing data.
Each record represents one exposure event.

Key Variables:
- PARAM/PARAMCD: Exposure parameter (e.g., "DOSE", "DURATION")
- AVAL: Analysis value
- ASTDT/AENDT: Start/end dates
- ADURN: Duration (numeric, in days)
"""
struct ADEXRecord
    # Identifiers
    studyid::String
    usubjid::String
    exseq::Int

    # Parameter
    paramcd::String   # Parameter code
    param::String     # Full parameter name
    paramn::Union{Nothing, Int}
    parcat1::String   # Category

    # Analysis values
    aval::Union{Nothing, Float64}     # Value (e.g., dose amount)
    avalu::String                      # Units
    avalc::String                      # Character value

    # Dates and duration
    astdt::Union{Nothing, Int}        # Start date (SAS date)
    astdtc::String                     # Start date (ISO 8601)
    aendt::Union{Nothing, Int}        # End date (SAS date)
    aendtc::String                     # End date (ISO 8601)
    adurn::Union{Nothing, Float64}    # Duration (days)
    aduru::String                      # Duration units

    # Treatment
    trtp::String
    trta::String

    # Administration
    exdosfrm::String  # Dose form
    exroute::String   # Route

    # Flags
    anl01fl::String

    # Visit
    avisit::String
    avisitn::Union{Nothing, Int}

    function ADEXRecord(;
        studyid::String="",
        usubjid::String="",
        exseq::Int=0,
        paramcd::String="DOSE",
        param::String="Dose",
        paramn::Union{Nothing, Int}=nothing,
        parcat1::String="",
        aval::Union{Nothing, Float64}=nothing,
        avalu::String="mg",
        avalc::String="",
        astdt::Union{Nothing, Int}=nothing,
        astdtc::String="",
        aendt::Union{Nothing, Int}=nothing,
        aendtc::String="",
        adurn::Union{Nothing, Float64}=nothing,
        aduru::String="DAYS",
        trtp::String="",
        trta::String="",
        exdosfrm::String="",
        exroute::String="",
        anl01fl::String="Y",
        avisit::String="",
        avisitn::Union{Nothing, Int}=nothing
    )
        new(studyid, usubjid, exseq, paramcd, param, paramn, parcat1,
            aval, avalu, avalc,
            astdt, astdtc, aendt, aendtc, adurn, aduru,
            trtp, trta, exdosfrm, exroute, anl01fl, avisit, avisitn)
    end
end

# ============================================================================
# Combined ADaM Dataset
# ============================================================================

"""
    ADaMDataset

Combined ADaM dataset containing all PK-relevant domains.

Fields:
- adsl: Subject-level analysis dataset
- adpc: PK concentration analysis dataset
- adpp: PK parameters analysis dataset
- adex: Exposure analysis dataset
- study_id: Study identifier
"""
struct ADaMDataset
    adsl::Vector{ADSLRecord}
    adpc::Vector{ADPCRecord}
    adpp::Vector{ADPPRecord}
    adex::Vector{ADEXRecord}
    study_id::String

    function ADaMDataset(;
        adsl::Vector{ADSLRecord}=ADSLRecord[],
        adpc::Vector{ADPCRecord}=ADPCRecord[],
        adpp::Vector{ADPPRecord}=ADPPRecord[],
        adex::Vector{ADEXRecord}=ADEXRecord[],
        study_id::String=""
    )
        new(adsl, adpc, adpp, adex, study_id)
    end
end

"""
    ADaMPKDataset

Simplified ADaM dataset focused on PK analysis.
Contains only ADSL (subjects), ADPC (concentrations), and ADEX (dosing).
"""
struct ADaMPKDataset
    adsl::Vector{ADSLRecord}
    adpc::Vector{ADPCRecord}
    adex::Vector{ADEXRecord}
    study_id::String

    function ADaMPKDataset(;
        adsl::Vector{ADSLRecord}=ADSLRecord[],
        adpc::Vector{ADPCRecord}=ADPCRecord[],
        adex::Vector{ADEXRecord}=ADEXRecord[],
        study_id::String=""
    )
        new(adsl, adpc, adex, study_id)
    end
end

# ============================================================================
# ADaM Standard Parameter Codes
# ============================================================================

"""
Standard PARAMCD values for ADPC (concentrations).
"""
const ADPC_PARAMCD = Dict(
    "concentration" => "CONC",
    "concentration_log" => "LCONC",
    "concentration_pred" => "PREDCONC",
    "concentration_ipred" => "IPREDCONC"
)

"""
Standard PARAMCD values for ADPP (PK parameters).
"""
const ADPP_PARAMCD = Dict(
    # Primary NCA parameters
    "auc_last" => "AUCLST",
    "auc_inf" => "AUCINF",
    "auc_0_t" => "AUC0T",
    "cmax" => "CMAX",
    "tmax" => "TMAX",
    "half_life" => "THALF",
    "lambda_z" => "LAMZ",
    "cl" => "CL",
    "vss" => "VSS",
    "vz" => "VZ",
    "mrt" => "MRT",

    # Secondary parameters
    "auc_0_12" => "AUC012",
    "auc_0_24" => "AUC024",
    "auc_extrap_pct" => "AUCPEXT",
    "cmin" => "CMIN",
    "cavg" => "CAVG",
    "ctrough" => "CTROUGH",
    "accumulation_ratio" => "RACC",

    # Compartmental parameters
    "ka" => "KA",
    "v1" => "V1",
    "v2" => "V2",
    "q" => "Q",
    "ke" => "KE"
)

export ADPC_PARAMCD, ADPP_PARAMCD

# ============================================================================
# Accessor Functions
# ============================================================================

"""Get unique subject IDs from ADaM dataset."""
function adam_subject_ids(data::ADaMDataset)::Vector{String}
    return unique([r.usubjid for r in data.adsl])
end

"""Get number of subjects in ADaM dataset."""
function adam_n_subjects(data::ADaMDataset)::Int
    return length(data.adsl)
end

"""Get all concentration records for a subject."""
function get_subject_adpc(data::ADaMDataset, usubjid::String)::Vector{ADPCRecord}
    return filter(r -> r.usubjid == usubjid, data.adpc)
end

"""Get all exposure records for a subject."""
function get_subject_adex(data::ADaMDataset, usubjid::String)::Vector{ADEXRecord}
    return filter(r -> r.usubjid == usubjid, data.adex)
end

"""Get ADSL record for a subject."""
function get_subject_adsl(data::ADaMDataset, usubjid::String)::Union{Nothing, ADSLRecord}
    idx = findfirst(r -> r.usubjid == usubjid, data.adsl)
    return idx === nothing ? nothing : data.adsl[idx]
end

export adam_subject_ids, adam_n_subjects, get_subject_adpc, get_subject_adex, get_subject_adsl
