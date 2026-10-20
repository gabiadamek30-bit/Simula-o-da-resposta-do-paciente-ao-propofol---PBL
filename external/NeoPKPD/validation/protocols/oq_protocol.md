# Operational Qualification (OQ) Protocol

**Document ID:** NEOPKPD-OQ-001
**Version:** 1.0
**Effective Date:** 2026-01-10
**Regulatory Reference:** FDA 21 CFR Part 11, GAMP 5

---

## 1. Purpose

This Operational Qualification (OQ) protocol verifies that NeoPKPD operates correctly within specified operating ranges and produces expected results. OQ demonstrates that the software functions as intended per the functional specifications.

## 2. Scope

This protocol covers:
- Core PK/PD simulation functionality
- Global sensitivity analysis (Sobol', Morris)
- Artifact serialization and reproducibility
- Compliance features (audit trail, integrity, environment capture)
- Golden artifact validation

## 3. Responsibilities

| Role | Responsibility |
|------|----------------|
| QA Lead | Protocol approval and final sign-off |
| Validation Engineer | OQ execution and documentation |
| Subject Matter Expert | Technical review of results |

## 4. Prerequisites

- IQ Protocol completed successfully
- Golden artifacts available in `validation/golden/`
- Test data sets available

## 5. Operational Verification Tests

### OQ-001: One-Compartment IV Bolus Simulation

**Objective:** Verify basic PK simulation produces expected results.

**Procedure:**
```julia
using NeoPKPD
spec = ModelSpec(
    OneCompIVBolus(),
    "test",
    OneCompIVBolusParams(1.0, 10.0),  # CL=1, V=10
    [DoseEvent(0.0, 100.0)]
)
grid = SimGrid(0.0, 24.0, collect(0.0:1.0:24.0))
solver = SolverSpec(:Tsit5, 1e-6, 1e-8, 10000)
result = simulate(spec, grid, solver)
```

**Expected Result:**
- Initial concentration C(0) = 100/10 = 10 mg/L
- Half-life = 0.693 * V / CL = 6.93 hours
- Exponential decay pattern

**Acceptance Criteria:** PASS if C(0) ≈ 10.0 ± 0.01

---

### OQ-002: Golden Artifact Reproducibility

**Objective:** Verify simulations can be exactly reproduced from golden artifacts.

**Procedure:**
```julia
using NeoPKPD, JSON

# Load golden artifact
artifact = JSON.parsefile("validation/golden/onecomp_iv_bolus_1.json")

# Replay and compare
original_result = artifact["result"]["concentration"]
replayed = replay_execution(artifact)
new_result = replayed["concentration"]
```

**Expected Result:** All concentration values match exactly

**Acceptance Criteria:** PASS if max absolute difference < 1e-10

---

### OQ-003: Sobol' Sensitivity Analysis

**Objective:** Verify GSA produces valid sensitivity indices.

**Procedure:**
```julia
using NeoPKPD
spec = ModelSpec(...)
bounds = ParameterBounds(Dict(:CL => (0.5, 2.0), :V => (5.0, 20.0)))
gsa_spec = GlobalSensitivitySpec(
    SobolMethod(base_sample_size=64, bootstrap_samples=0),
    bounds
)
result = run_sobol_sensitivity(spec, grid, solver, gsa_spec)
```

**Expected Result:**
- 0 ≤ Si ≤ 1 for all parameters
- 0 ≤ STi ≤ 1 for all parameters
- STi ≥ Si for all parameters

**Acceptance Criteria:** PASS if all indices within valid ranges

---

### OQ-004: Morris Elementary Effects

**Objective:** Verify Morris screening method works correctly.

**Procedure:**
```julia
using NeoPKPD
gsa_spec = GlobalSensitivitySpec(
    MorrisMethod(n_trajectories=10),
    bounds
)
result = run_morris_sensitivity(spec, grid, solver, gsa_spec)
```

**Expected Result:**
- μ* ≥ 0 for all parameters
- σ ≥ 0 for all parameters
- Elementary effects computed for each parameter

**Acceptance Criteria:** PASS if indices computed without error

---

### OQ-005: Artifact Serialization Round-Trip

**Objective:** Verify artifacts can be serialized and deserialized without data loss.

**Procedure:**
```julia
using NeoPKPD, JSON

# Serialize
artifact = serialize_execution(model_spec=spec, grid=grid, solver=solver, result=result)
json_str = JSON.json(artifact)

# Deserialize
artifact2 = JSON.parse(json_str)
parsed = deserialize_execution(artifact2)

# Replay
result2 = replay_execution(artifact2)
```

**Expected Result:** Replayed result matches original

**Acceptance Criteria:** PASS if results identical

---

### OQ-006: Compliance Metadata Generation

**Objective:** Verify compliance metadata is correctly added to artifacts.

**Procedure:**
```julia
using NeoPKPD

artifact = serialize_execution(...)
add_compliance_metadata!(artifact)

meta = artifact["compliance_metadata"]
```

**Expected Result:**
- `audit_record` contains valid UUID and timestamp
- `integrity` contains valid SHA-256 hashes
- `environment` contains valid system information

**Acceptance Criteria:** PASS if all compliance fields present and valid

---

### OQ-007: Integrity Verification

**Objective:** Verify artifact integrity can be validated.

**Procedure:**
```julia
using NeoPKPD

artifact = serialize_execution(...)
add_compliance_metadata!(artifact)

# Verify unmodified artifact
result1 = verify_artifact_integrity(artifact)

# Modify and verify
artifact["result"]["concentration"][1] = 999.0
result2 = verify_artifact_integrity(artifact)
```

**Expected Result:**
- `result1.is_valid == true`
- `result2.is_valid == false`

**Acceptance Criteria:** PASS if integrity correctly detected

---

### OQ-008: Audit Trail Verification

**Objective:** Verify audit records are self-verifying.

**Procedure:**
```julia
using NeoPKPD

record = create_audit_record()
is_valid = verify_audit_record(record)
```

**Expected Result:** `is_valid == true`

**Acceptance Criteria:** PASS if audit record verifies

---

### OQ-009: Reproducibility with Seed

**Objective:** Verify simulations with same seed produce identical results.

**Procedure:**
```julia
using NeoPKPD

gsa_spec1 = GlobalSensitivitySpec(..., seed=UInt64(42))
gsa_spec2 = GlobalSensitivitySpec(..., seed=UInt64(42))

result1 = run_sobol_sensitivity(spec, grid, solver, gsa_spec1)
result2 = run_sobol_sensitivity(spec, grid, solver, gsa_spec2)
```

**Expected Result:** Results are bitwise identical

**Acceptance Criteria:** PASS if all indices match exactly

---

### OQ-010: Golden Artifact Integrity

**Objective:** Verify all golden artifacts pass integrity checks.

**Procedure:**
```julia
using NeoPKPD, JSON

for file in readdir("validation/golden")
    artifact = JSON.parsefile(joinpath("validation/golden", file))
    if has_compliance_metadata(artifact)
        result = verify_artifact_integrity(artifact)
        @assert result.is_valid "Failed: $file"
    end
end
```

**Expected Result:** All artifacts pass integrity verification

**Acceptance Criteria:** PASS if all golden artifacts valid

---

## 6. Results Summary

| Test ID | Description | Result | Tester | Date |
|---------|-------------|--------|--------|------|
| OQ-001 | IV Bolus Simulation | ☐ PASS ☐ FAIL | | |
| OQ-002 | Golden Reproducibility | ☐ PASS ☐ FAIL | | |
| OQ-003 | Sobol' Analysis | ☐ PASS ☐ FAIL | | |
| OQ-004 | Morris Analysis | ☐ PASS ☐ FAIL | | |
| OQ-005 | Serialization Round-Trip | ☐ PASS ☐ FAIL | | |
| OQ-006 | Compliance Metadata | ☐ PASS ☐ FAIL | | |
| OQ-007 | Integrity Verification | ☐ PASS ☐ FAIL | | |
| OQ-008 | Audit Trail | ☐ PASS ☐ FAIL | | |
| OQ-009 | Seed Reproducibility | ☐ PASS ☐ FAIL | | |
| OQ-010 | Golden Integrity | ☐ PASS ☐ FAIL | | |

## 7. Deviations and Corrective Actions

| Deviation # | Test ID | Description | Corrective Action | Status |
|-------------|---------|-------------|-------------------|--------|
| | | | | |

## 8. Conclusion

☐ All OQ tests PASSED - Operations are qualified
☐ Deviations documented and resolved
☐ OQ Protocol execution COMPLETE

---

## 9. Approval Signatures

| Role | Name | Signature | Date |
|------|------|-----------|------|
| Executed By | | | |
| Reviewed By | | | |
| Approved By | | | |

---

*Document generated by NeoPKPD Compliance Module*
