# Contributing to NeoPKPD

Thank you for your interest in contributing to NeoPKPD! This document provides guidelines and information for contributors.

## Code of Conduct

Please be respectful and constructive in all interactions. We welcome contributors of all backgrounds and experience levels.

## Core Principles

NeoPKPD maintains strict standards for scientific computing:

1. **Deterministic Execution**: Same inputs must always produce identical outputs
2. **Explicit Configuration**: No hidden defaults or implicit behavior
3. **Versioned Semantics**: All numerical behavior changes require version bumps
4. **Complete Traceability**: Every simulation can be fully reproduced from artifacts

## Getting Started

### Prerequisites

- Julia 1.9+
- Python 3.9+ (for Python bindings)
- Git

### Development Setup

```bash
# Clone the repository
git clone https://github.com/neopkpd/neopkpd.git
cd neopkpd

# Install Julia dependencies
julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project=packages/core -e 'using Pkg; Pkg.test()'

# (Optional) Set up Python development
cd python
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
pytest tests/
```

## Repository Structure

```
neopkpd/
├── packages/core/    # Core Julia package
├── cli/NeoPKPDCLI/      # Command-line interface
├── python/               # Python bindings
├── docs/                 # Documentation (MkDocs)
├── validation/           # Golden artifacts and validation
├── scripts/              # Development scripts
└── bin/                  # CLI entry point
```

## Making Changes

### Branch Naming

- `feature/description` - New features
- `fix/description` - Bug fixes
- `docs/description` - Documentation changes
- `refactor/description` - Code refactoring

### Code Style

**Julia:**
- Use JuliaFormatter with project settings
- Run: `julia scripts/format.jl`

**Python:**
- Follow PEP 8
- Use type hints where appropriate

### Commit Messages

Follow conventional commits:

```
type(scope): description

[optional body]

[optional footer]
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

## Pull Request Process

### Requirements

Every PR must include:

1. **Tests** covering new behavior
2. **Documentation** for new features
3. **Golden artifact updates** if numerical outputs change
4. **Changelog entry** for user-facing changes

### Checklist

- [ ] Tests pass: `julia --project=packages/core -e 'using Pkg; Pkg.test()'`
- [ ] Golden validation passes: `./bin/neopkpd validate-golden`
- [ ] Documentation builds: `mkdocs build --strict`
- [ ] Python tests pass (if applicable): `cd python && pytest tests/`
- [ ] Code is formatted

### Review Process

1. Open PR with clear description
2. CI checks must pass
3. At least one maintainer review required
4. Address all feedback
5. Squash merge to main

## Semantic Versioning

NeoPKPD uses three semantic versions:

| Version | Scope | Bump When |
|---------|-------|-----------|
| Event Semantics | Dose handling | Dose normalization changes |
| Solver Semantics | ODE solving | Solver behavior changes |
| Artifact Schema | JSON format | Artifact structure changes |

### Bumping Versions

If your change affects numerical output:

1. Increment appropriate version in `packages/core/src/engine/`
2. Regenerate golden artifacts: `julia validation/scripts/generate_golden_artifacts.jl`
3. Document the change in CHANGELOG.md

## Adding New Features

### New PK Model

1. Create `packages/core/src/models/new_model.jl`
2. Implement required interface (see `pk_interface.jl`)
3. Add serialization support
4. Add tests
5. Create golden artifact
6. Document in `docs/models.md`

### New PD Model

1. Create `packages/core/src/pd/new_pd.jl`
2. Implement required interface
3. Add serialization support
4. Add tests
5. Document in `docs/models.md`

## Reporting Issues

### Bug Reports

Include:
- NeoPKPD version
- Julia/Python version
- Minimal reproduction code
- Expected vs actual behavior
- Error messages (full traceback)

### Feature Requests

Include:
- Use case description
- Proposed API (if applicable)
- Scientific rationale

## Questions?

- Open a GitHub Discussion for questions
- Open an Issue for bugs or feature requests
- Check existing documentation first

## License

By contributing, you agree that your contributions will be licensed under the project's MIT License.
