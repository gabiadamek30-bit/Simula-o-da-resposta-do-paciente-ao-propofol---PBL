"""
Pytest configuration for NeoPKPD Python tests.

This module provides shared fixtures and configuration for all tests,
including proper Julia initialization with signal handling to prevent
segmentation faults in CI environments.
"""

import os
import pytest

# Configure juliacall signal handling BEFORE any Julia imports
# This must happen before juliacall is imported
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")


@pytest.fixture(scope="session", autouse=True)
def julia_session():
    """
    Session-scoped fixture that initializes Julia once for all tests.

    This fixture:
    1. Sets up proper signal handling for juliacall
    2. Initializes the Julia runtime before any tests run
    3. Ensures Julia is only initialized once per test session

    The autouse=True ensures this runs before any test, preventing
    segmentation faults that can occur when Julia is initialized
    mid-test-run in CI environments.
    """
    import neopkpd

    # Initialize Julia with proper signal handling
    try:
        neopkpd.init_julia()
    except Exception as e:
        pytest.skip(f"Failed to initialize Julia: {e}")

    yield

    # Cleanup (if needed in future)


@pytest.fixture(scope="module")
def init():
    """
    Module-scoped fixture for tests that need Julia initialized.

    This is a backward-compatible fixture that tests can use to ensure
    Julia is available. The actual initialization happens in julia_session.
    """
    import neopkpd

    # Julia should already be initialized by julia_session
    # This just ensures it's available
    try:
        neopkpd._core._require_julia()
    except Exception as e:
        pytest.skip(f"Julia not available: {e}")

    yield


@pytest.fixture(scope="module")
def julia_initialized():
    """
    Module-scoped fixture for tests that need Julia initialized.

    Alias for 'init' fixture for backward compatibility with existing tests.
    The actual initialization happens in julia_session (autouse=True).
    """
    import neopkpd

    # Julia should already be initialized by julia_session
    try:
        neopkpd._core._require_julia()
    except Exception as e:
        pytest.skip(f"Julia not available: {e}")

    return True


@pytest.fixture(scope="function")
def julia():
    """
    Function-scoped fixture that provides access to Julia main module.

    Usage:
        def test_something(julia):
            result = julia.NeoPKPD.some_function()
    """
    import neopkpd

    jl = neopkpd._core._require_julia()
    if jl is None:
        pytest.skip("Julia not initialized")

    return jl


# Skip markers for tests that require specific features
def pytest_configure(config):
    """Configure custom pytest markers."""
    config.addinivalue_line(
        "markers", "julia: mark test as requiring Julia runtime"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to handle Julia-dependent tests."""
    # Check if Julia is available
    julia_available = True
    try:
        import neopkpd
        neopkpd.init_julia()
    except Exception:
        julia_available = False

    if not julia_available:
        skip_julia = pytest.mark.skip(reason="Julia runtime not available")
        for item in items:
            # Skip tests that use the 'init' or 'julia' fixtures
            if "init" in item.fixturenames or "julia" in item.fixturenames:
                item.add_marker(skip_julia)
