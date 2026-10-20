# Release procedure

Releases are intentional and conservative.

## Steps

1. Ensure main branch is green.
2. Decide release type (PATCH, MINOR, MAJOR).
3. Update:
   - VERSION
   - packages/core/src/NeoPKPD.jl
   - CHANGELOG.md
4. If MAJOR:
   - Bump semantics version constants
   - Regenerate golden artifacts
5. Commit with message: "Release X.Y.Z"
6. Tag:
   git tag -a vX.Y.Z -m "NeoPKPD X.Y.Z"
7. Push tag:
   git push origin vX.Y.Z
8. Create GitHub release from tag with CHANGELOG excerpt.

## Artifact compatibility

Artifacts are guaranteed replayable within the same MAJOR version.
Cross-MAJOR replay requires explicit migration support.
