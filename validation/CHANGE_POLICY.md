# Change policy

NeoPKPD treats numerical behavior as versioned scientific infrastructure.

## Golden artifacts

Golden artifacts in validation/golden define expected numerical behavior.

Any change in:
- states
- observations
- time grid

is considered a breaking numerical change.

## Required action on change

If golden artifacts change, at least one of the following must be bumped:

- EVENT_SEMANTICS_VERSION
- SOLVER_SEMANTICS_VERSION
- ARTIFACT_SCHEMA_VERSION

The bump must be intentional and documented in the commit message.

## Non breaking changes

The following do not require a semantics bump:
- Documentation changes
- Refactoring that preserves golden equality
- Performance improvements with identical outputs

## Rationale

This policy ensures:
- Long-term reproducibility
- Trustworthy regression detection
- Explicit scientific intent
