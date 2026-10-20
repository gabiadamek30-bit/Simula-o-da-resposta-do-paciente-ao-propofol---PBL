# Release Process

This document describes how to release new versions of NeoPKPD (Julia) and neopkpd (Python).

## Version Numbers

Both packages use the same version number and should be released together:
- Julia: `packages/core/Project.toml`
- Python: `packages/python/pyproject.toml`

## Pre-Release Checklist

Before releasing, ensure:

- [ ] All tests pass: `julia --project=packages/core -e 'using Pkg; Pkg.test()'`
- [ ] Python tests pass: `cd packages/python && pytest tests/`
- [ ] CHANGELOG.md is updated with all changes
- [ ] Version numbers are consistent across all files
- [ ] Documentation builds: `mkdocs build --strict`
- [ ] Golden artifact validation passes

## Release Steps

### 1. Update Version Numbers

Update the version in both packages:

```bash
# packages/core/Project.toml
version = "X.Y.Z"

# packages/python/pyproject.toml
version = "X.Y.Z"
```

### 2. Update CHANGELOG

Add a new section to `CHANGELOG.md` with the release date:

```markdown
## [X.Y.Z] - YYYY-MM-DD

### Added
- New features...

### Changed
- Changes...

### Fixed
- Bug fixes...
```

### 3. Commit Changes

```bash
git add -A
git commit -m "Release v X.Y.Z"
```

### 4. Create and Push Tag

```bash
git tag -a vX.Y.Z -m "Release vX.Y.Z"
git push origin main
git push origin vX.Y.Z
```

This will trigger the release workflow which will:
- Run all tests
- Build the Python package
- Publish to PyPI
- Create a GitHub release

### 5. Register Julia Package (First Time)

For the first release, you need to register with the Julia General Registry:

1. Install Registrator:
   ```julia
   using Pkg
   Pkg.add("Registrator")
   ```

2. Comment on the release commit or GitHub release:
   ```
   @JuliaRegistrator register subdir=packages/core
   ```

3. Wait for the registry PR to be merged (may take 3 days for new packages)

For subsequent releases, just comment `@JuliaRegistrator register subdir=packages/core` on the new tag.

## Manual Release (if needed)

### Python Package to PyPI

```bash
cd packages/python
pip install build twine

# Build
python -m build

# Upload to TestPyPI first
twine upload --repository testpypi dist/*

# Upload to PyPI
twine upload dist/*
```

### Dry Run

You can test the release workflow without publishing:

1. Go to Actions → Release workflow
2. Click "Run workflow"
3. Enter the version number
4. Check "Dry run"
5. Click "Run workflow"

This will run all tests and build packages but publish to TestPyPI instead of PyPI.

## PyPI Configuration

To publish to PyPI, you need to configure:

1. Create a PyPI account at https://pypi.org
2. Create an API token with upload permissions
3. Add the token as a GitHub secret named `PYPI_API_TOKEN`

Or use Trusted Publishing (recommended):
1. Go to PyPI → Your projects → neopkpd → Publishing
2. Add a new publisher with:
   - Owner: `shramish2057`
   - Repository: `openpkpd`
   - Workflow: `release.yml`
   - Environment: `pypi`

## Versioning Policy

We follow [Semantic Versioning](https://semver.org/):

- **MAJOR** (X.0.0): Breaking API changes
- **MINOR** (0.X.0): New features, backwards compatible
- **PATCH** (0.0.X): Bug fixes, backwards compatible

Additionally, we version numerical semantics independently:
- Event Semantics: Changes to dose/event handling
- Solver Semantics: Changes to ODE solver behavior
- Artifact Schema: Changes to JSON format

Any change to numerical output requires incrementing the appropriate semantics version.

## Post-Release

After a successful release:

1. Announce on relevant channels
2. Update documentation if needed
3. Close related GitHub issues/milestones
4. Start the next development cycle

## Troubleshooting

### PyPI upload fails
- Check that the version doesn't already exist on PyPI
- Verify API token permissions
- Try uploading to TestPyPI first

### Julia registration fails
- Ensure `Project.toml` has valid UUID
- Check that version follows semver
- Verify package name doesn't conflict with existing packages

### Tests fail on release
- Fix the issues and create a new patch release
- Never force-push tags that have been pushed to remote
