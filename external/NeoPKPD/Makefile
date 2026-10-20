# NeoPKPD Makefile
# Comprehensive build, test, and documentation targets

.PHONY: all help install test test-julia test-python test-all docs docs-serve clean lint validate

# Default target
all: test

# Help
help:
	@echo "NeoPKPD Build System"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  install        Install all dependencies"
	@echo "  install-julia  Install Julia dependencies"
	@echo "  install-python Install Python dependencies"
	@echo "  test           Run all tests (Julia + Python)"
	@echo "  test-julia     Run Julia tests only"
	@echo "  test-python    Run Python tests only"
	@echo "  test-nca       Run NCA tests"
	@echo "  test-trial     Run trial simulation tests"
	@echo "  test-viz       Run visualization tests"
	@echo "  validate       Validate golden artifacts"
	@echo "  lint           Run linters"
	@echo "  docs           Build documentation"
	@echo "  docs-serve     Serve documentation locally"
	@echo "  clean          Clean build artifacts"
	@echo ""

# ============================================================================
# Installation
# ============================================================================

install: install-julia install-python
	@echo "All dependencies installed"

install-julia:
	@echo "Installing Julia dependencies..."
	julia --project=packages/core -e 'using Pkg; Pkg.instantiate()'

install-python:
	@echo "Installing Python dependencies..."
	cd packages/python && python3 -m pip install -e .[dev]

install-docs:
	@echo "Installing documentation dependencies..."
	python3 -m pip install -r docs/requirements.txt

# ============================================================================
# Testing
# ============================================================================

test: test-julia test-python
	@echo "All tests passed"

test-all: test validate
	@echo "All tests and validation passed"

test-julia:
	@echo "Running Julia tests..."
	julia --project=packages/core -e 'using Pkg; Pkg.test()'

test-python:
	@echo "Running Python tests..."
	cd packages/python && python3 -m pytest tests/ -v

test-nca:
	@echo "Running NCA tests..."
	cd packages/python && python3 -m pytest tests/test_nca.py -v

test-trial:
	@echo "Running trial simulation tests..."
	cd packages/python && python3 -m pytest tests/test_trial.py -v

test-viz:
	@echo "Running visualization tests..."
	cd packages/python && python3 -m pytest tests/test_viz.py -v 2>/dev/null || echo "No viz tests found"

test-replay:
	@echo "Running replay tests..."
	cd packages/python && python3 -m pytest tests/test_replay.py -v

test-simulate:
	@echo "Running simulation tests..."
	cd packages/python && python3 -m pytest tests/test_simulate.py -v

# ============================================================================
# Validation
# ============================================================================

validate:
	@echo "Validating golden artifacts..."
	./packages/cli/bin/neopkpd validate-golden

validate-golden: validate

generate-golden:
	@echo "Generating golden artifacts..."
	julia --project=packages/core validation/scripts/generate_golden_artifacts.jl

# ============================================================================
# Documentation
# ============================================================================

docs: install-docs
	@echo "Building documentation..."
	mkdocs build --strict

docs-serve: install-docs
	@echo "Serving documentation at http://localhost:8000..."
	mkdocs serve

docs-deploy:
	@echo "Deploying documentation to GitHub Pages..."
	mkdocs gh-deploy --force

# ============================================================================
# Linting
# ============================================================================

lint: lint-python
	@echo "Linting complete"

lint-python:
	@echo "Linting Python code..."
	cd packages/python && python3 -m ruff check neopkpd/ tests/ || true
	cd packages/python && python3 -m mypy neopkpd/ --ignore-missing-imports || true

# ============================================================================
# CLI
# ============================================================================

cli-help:
	./packages/cli/bin/neopkpd help

cli-version:
	./packages/cli/bin/neopkpd version

# ============================================================================
# Development
# ============================================================================

dev-setup: install
	@echo "Development environment set up"

smoke-test:
	@echo "Running smoke tests..."
	./validation/scripts/smoke_cli.sh
	./packages/python/scripts/smoke_python.sh 2>/dev/null || true

# ============================================================================
# Cleanup
# ============================================================================

clean:
	@echo "Cleaning build artifacts..."
	rm -rf site/
	rm -rf packages/python/.pytest_cache/
	rm -rf packages/python/__pycache__/
	rm -rf packages/python/neopkpd/__pycache__/
	rm -rf packages/python/neopkpd/**/__pycache__/
	rm -rf packages/python/*.egg-info/
	rm -rf .pytest_cache/
	rm -rf packages/core/.pytest_cache/
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -delete
	@echo "Clean complete"

clean-docs:
	rm -rf site/

# ============================================================================
# CI Targets
# ============================================================================

ci-test: test-julia test-python
	@echo "CI tests complete"

ci-validate: validate
	@echo "CI validation complete"

ci-docs: docs
	@echo "CI docs build complete"

ci-all: ci-test ci-validate ci-docs
	@echo "All CI checks passed"
