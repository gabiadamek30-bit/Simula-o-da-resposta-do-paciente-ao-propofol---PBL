# Performance Qualification (PQ) Protocol

**Document ID:** NEOPKPD-PQ-001
**Version:** 1.0
**Effective Date:** 2026-01-10
**Regulatory Reference:** FDA 21 CFR Part 11, GAMP 5

---

## 1. Purpose

This Performance Qualification (PQ) protocol verifies that NeoPKPD performs consistently and reliably under actual operating conditions. PQ demonstrates that the system meets acceptance criteria during real-world use cases.

## 2. Scope

This protocol covers:
- End-to-end simulation workflows
- Population PK/PD simulations
- Clinical trial simulations
- Performance under load
- Cross-platform consistency

## 3. Responsibilities

| Role | Responsibility |
|------|----------------|
| QA Lead | Protocol approval and final sign-off |
| Pharmacometrician | Scenario design and result validation |
| Validation Engineer | PQ execution and documentation |

## 4. Prerequisites

- IQ Protocol completed successfully
- OQ Protocol completed successfully
- Representative test scenarios available
- Baseline performance metrics established

## 5. Performance Verification Tests

### PQ-001: End-to-End Workflow - Drug Development

**Objective:** Verify complete simulation workflow for a typical drug development scenario.

**Scenario:**
- One-compartment IV bolus model with IIV
- Population of 100 virtual subjects
- Multiple dosing regimen
- PK endpoint: AUC, Cmax, Tmax

**Procedure:**
```julia
using NeoPKPD

# Define population
iiv = IIV([
    IIVParameter(:CL, LogNormal(), 0.3),
    IIVParameter(:V, LogNormal(), 0.25)
])

pop_spec = PopulationSpec(
    OneCompIVBolus(),
    "phase1_sim",
    OneCompIVBolusParams(1.0, 10.0),
    [DoseEvent(0.0, 100.0), DoseEvent(24.0, 100.0)],
    100,
    iiv
)

grid = SimGrid(0.0, 72.0, collect(0.0:0.5:72.0))
solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)

result = simulate_population(pop_spec, grid, solver)
```

**Expected Result:**
- All 100 subjects simulated
- Variability in PK profiles reflects IIV
- AUC and Cmax vary across subjects

**Acceptance Criteria:** PASS if simulation completes and results are pharmacologically plausible

---

### PQ-002: Sensitivity Analysis Workflow

**Objective:** Verify GSA workflow produces actionable insights.

**Scenario:**
- Two-compartment model with 5 parameters
- Sobol' analysis with bootstrap CI
- Morris screening for comparison

**Procedure:**
```julia
using NeoPKPD

# Run Sobol'
bounds = ParameterBounds(Dict(
    :CL => (0.5, 2.0),
    :V1 => (5.0, 20.0),
    :Q => (0.1, 1.0),
    :V2 => (10.0, 50.0),
    :ka => (0.5, 2.0)
))

sobol_result = run_sobol_sensitivity(spec, grid, solver,
    GlobalSensitivitySpec(
        SobolMethod(base_sample_size=128, bootstrap_samples=100),
        bounds
    )
)

morris_result = run_morris_sensitivity(spec, grid, solver,
    GlobalSensitivitySpec(
        MorrisMethod(n_trajectories=20),
        bounds
    )
)

# Compare rankings
sobol_ranking = rank_parameters(sobol_result; by=:STi)
morris_ranking = rank_parameters(morris_result; by=:mu_star)
```

**Expected Result:**
- Both methods identify similar important parameters
- Confidence intervals provided for Sobol' indices
- Results aid in parameter prioritization

**Acceptance Criteria:** PASS if top 2 parameters consistent between methods

---

### PQ-003: Clinical Trial Simulation

**Objective:** Verify trial simulation produces regulatory-grade outputs.

**Scenario:**
- Parallel group design (placebo, low dose, high dose)
- 50 subjects per arm
- 12-week treatment period
- Power analysis for primary endpoint

**Expected Result:**
- Complete trial simulation
- Summary statistics per arm
- Confidence intervals for treatment effect

**Acceptance Criteria:** PASS if trial completes with valid statistics

---

### PQ-004: Artifact Audit Trail Workflow

**Objective:** Verify audit trail captures complete execution history.

**Scenario:**
1. Create initial simulation
2. Serialize with compliance metadata
3. Modify parameters and re-run
4. Serialize as replay
5. Verify audit chain

**Procedure:**
```julia
using NeoPKPD

# Initial execution
result1 = simulate(spec, grid, solver)
artifact1 = serialize_execution(...)
add_compliance_metadata!(artifact1, action=AUDIT_CREATE)

# Replay/modify
artifact2 = serialize_execution(...)
add_compliance_metadata!(artifact2,
    action=AUDIT_REPLAY,
    previous_execution_id=get_execution_id(artifact1)
)

# Verify chain
@assert artifact2["compliance_metadata"]["audit_record"]["previous_execution_id"] ==
        artifact1["compliance_metadata"]["audit_record"]["execution_id"]
```

**Expected Result:**
- Audit chain maintained between artifacts
- All timestamps in UTC
- All execution IDs unique

**Acceptance Criteria:** PASS if audit chain verified

---

### PQ-005: Cross-Platform Reproducibility

**Objective:** Verify results are reproducible across platforms.

**Scenario:**
- Run same simulation on macOS, Linux, Windows
- Same random seed
- Compare results

**Expected Result:**
- Numerical results identical (or within tolerance)
- Environment differences captured in metadata

**Acceptance Criteria:** PASS if max relative difference < 1e-10

---

### PQ-006: Performance Benchmarks

**Objective:** Verify system meets performance requirements.

**Test Cases:**
| Scenario | Expected Time | Tolerance |
|----------|---------------|-----------|
| Single subject, 1000 time points | < 1 sec | +50% |
| Population 100 subjects | < 10 sec | +50% |
| Sobol' N=1000 | < 60 sec | +100% |
| Morris r=50 | < 30 sec | +100% |

**Acceptance Criteria:** PASS if all benchmarks within tolerance

---

### PQ-007: Error Handling and Recovery

**Objective:** Verify graceful handling of edge cases.

**Test Cases:**
1. Invalid parameter values (negative CL)
2. Empty dose list
3. Zero time grid
4. Solver convergence failure

**Expected Result:**
- Clear error messages
- No data corruption
- Recoverable state

**Acceptance Criteria:** PASS if all errors handled gracefully

---

### PQ-008: Validation Report Generation

**Objective:** Verify automated validation reports are generated correctly.

**Procedure:**
```julia
using NeoPKPD

# Run IQ tests
include("validation/scripts/run_iq.jl")

# Run OQ tests
include("validation/scripts/run_oq.jl")

# Verify report generation
@assert isfile("validation/reports/iq_report.md")
@assert isfile("validation/reports/oq_report.md")
```

**Expected Result:**
- Reports generated in Markdown format
- All test results documented
- Environment captured in report

**Acceptance Criteria:** PASS if reports generated and contain all test results

---

## 6. Results Summary

| Test ID | Description | Result | Tester | Date |
|---------|-------------|--------|--------|------|
| PQ-001 | E2E Drug Development | ☐ PASS ☐ FAIL | | |
| PQ-002 | GSA Workflow | ☐ PASS ☐ FAIL | | |
| PQ-003 | Trial Simulation | ☐ PASS ☐ FAIL | | |
| PQ-004 | Audit Trail | ☐ PASS ☐ FAIL | | |
| PQ-005 | Cross-Platform | ☐ PASS ☐ FAIL ☐ N/A | | |
| PQ-006 | Performance | ☐ PASS ☐ FAIL | | |
| PQ-007 | Error Handling | ☐ PASS ☐ FAIL | | |
| PQ-008 | Report Generation | ☐ PASS ☐ FAIL | | |

## 7. Performance Metrics

| Metric | Baseline | Current | Status |
|--------|----------|---------|--------|
| Single subject simulation | | | |
| Population simulation (N=100) | | | |
| Sobol' analysis (N=1000) | | | |
| Morris screening (r=50) | | | |

## 8. Deviations and Corrective Actions

| Deviation # | Test ID | Description | Corrective Action | Status |
|-------------|---------|-------------|-------------------|--------|
| | | | | |

## 9. Conclusion

☐ All PQ tests PASSED - System performance qualified
☐ Deviations documented and resolved
☐ PQ Protocol execution COMPLETE
☐ System approved for production use

---

## 10. Approval Signatures

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Executed By | | | |
| SME Review | | | |
| QA Approval | | | |
| Final Approval | | | |

---

*Document generated by NeoPKPD Compliance Module*
