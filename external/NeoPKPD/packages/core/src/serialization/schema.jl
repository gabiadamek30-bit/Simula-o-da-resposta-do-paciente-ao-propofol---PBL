export ARTIFACT_SCHEMA_VERSION

# Schema version history:
# 1.0.0 - Initial release
# 1.1.0 - (planned) compliance_metadata block for FDA 21 CFR Part 11 alignment
#         - audit_record: execution ID, timestamps, user/host info
#         - integrity: SHA-256 hashes for content, inputs, outputs
#         - environment: complete environment snapshot
#         (Not yet activated - requires golden artifact regeneration)
const ARTIFACT_SCHEMA_VERSION = "1.0.0"
