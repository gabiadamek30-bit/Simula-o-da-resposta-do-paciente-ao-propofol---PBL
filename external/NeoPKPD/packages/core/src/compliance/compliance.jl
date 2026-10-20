# FDA 21 CFR Part 11 Compliance Module
# =====================================
#
# This module provides regulatory compliance features aligned with:
# - FDA 21 CFR Part 11: Electronic Records and Signatures
# - GAMP 5: Good Automated Manufacturing Practice
# - ALCOA+ Principles: Attributable, Legible, Contemporaneous, Original, Accurate,
#                      Complete, Consistent, Enduring, Available
#
# Features:
# - Audit trails with execution IDs and timestamps
# - Data integrity via SHA-256 checksums
# - Environment capture for reproducibility
# - Validation documentation (IQ/OQ/PQ)
#
# Usage:
# ```julia
# using NeoPKPD
#
# # Default: COMPLIANCE_STANDARD level
# result = simulate(spec, grid, solver)
# write_execution_json("output.json"; ...) # Includes compliance_metadata
#
# # Disable for development
# with_compliance(COMPLIANCE_DISABLED) do
#     result = simulate(spec, grid, solver)
# end
#
# # Verify artifact integrity
# verify_artifact_integrity(artifact)
# ```

using Pkg
using Dates

# Configuration and levels
include("config.jl")

# Environment capture
include("environment.jl")

# Integrity (SHA-256 hashing) - Phase 2
include("integrity.jl")

# Audit trail - Phase 2
include("audit_trail.jl")

# Digital signatures - Phase 1 Gold-Tier
include("signatures.jl")

# JSON Schema validation - Phase 1 Gold-Tier
include("schema_validation.jl")

# Schema migration framework - Phase 1 Gold-Tier
include("migration.jl")

# Validation reports - Phase 4
include("validation_report.jl")
