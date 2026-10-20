#!/usr/bin/env python3
"""
NeoPKPD Benchmark Suite - Python Bindings
==========================================

Benchmarks NeoPKPD Python interface (via juliacall).
Measures overhead of Python-Julia bridge.

Usage:
    cd packages/python
    python ../core/benchmarks/scripts/benchmark_neopkpd_python.py

Output:
    benchmarks/results/neopkpd_python_benchmarks.csv
"""

import os
import sys
import time
import statistics
from datetime import datetime
from pathlib import Path
import csv

# Set juliacall signal handling before any imports
os.environ.setdefault("PYTHON_JULIACALL_HANDLE_SIGNALS", "yes")

# =============================================================================
# Configuration
# =============================================================================

CONFIG = {
    "n_runs": 100,
    "n_warmup": 5,
    "random_seed": 12345,
    "output_dir": Path(__file__).parent.parent / "results",
}

# =============================================================================
# Setup
# =============================================================================

print("\n" + "=" * 70)
print("NeoPKPD Python Benchmark Suite")
print("=" * 70)
print(f"Python Version: {sys.version}")
print(f"Runs per benchmark: {CONFIG['n_runs']}")
print(f"Timestamp: {datetime.now()}")
print("=" * 70)
print("\nInitializing Julia runtime (this may take a moment)...")

# Initialize NeoPKPD
import neopkpd
neopkpd.init_julia()

print(f"NeoPKPD Version: {neopkpd.version()}")
print("Julia initialized successfully.\n")

# =============================================================================
# Benchmark Infrastructure
# =============================================================================

class BenchmarkResult:
    def __init__(self, category, name, model, n_subjects, times):
        self.category = category
        self.name = name
        self.model = model
        self.n_subjects = n_subjects
        self.times = times
        self.timestamp = datetime.now().isoformat()
        self.python_version = sys.version.split()[0]
        self.neopkpd_version = neopkpd.version()

    def summarize(self):
        times_ms = [t * 1000 for t in self.times]
        sorted_times = sorted(times_ms)
        n = len(times_ms)

        return {
            "category": self.category,
            "name": self.name,
            "model": self.model,
            "n_subjects": self.n_subjects,
            "n_runs": n,
            "mean_ms": statistics.mean(times_ms),
            "std_ms": statistics.stdev(times_ms) if n > 1 else 0,
            "median_ms": statistics.median(times_ms),
            "min_ms": min(times_ms),
            "max_ms": max(times_ms),
            "ci_lower_ms": sorted_times[int(n * 0.025)],
            "ci_upper_ms": sorted_times[int(n * 0.975)],
            "timestamp": self.timestamp,
            "python_version": self.python_version,
            "neopkpd_version": self.neopkpd_version,
        }


def run_benchmark(func, name, category, model, n_subjects=1):
    """Run a benchmark with warmup and timing."""
    print(f"  Running: {name} ({model})...")

    # Warmup
    for _ in range(CONFIG["n_warmup"]):
        func()

    # Timed runs
    times = []
    for _ in range(CONFIG["n_runs"]):
        start = time.perf_counter()
        func()
        end = time.perf_counter()
        times.append(end - start)

    result = BenchmarkResult(category, name, model, n_subjects, times)
    s = result.summarize()
    print(f"    Mean: {s['mean_ms']:.3f} ms (±{s['std_ms']:.3f}), Median: {s['median_ms']:.3f} ms")

    return result


# =============================================================================
# Benchmark Definitions
# =============================================================================

def benchmark_single_simulations():
    """Single simulation benchmarks."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Single Simulation Speed (Python)")
    print("=" * 60)

    results = []
    grid = {"tspan": (0.0, 24.0), "saveat": [i * 0.5 for i in range(49)]}
    solver = {"algorithm": "Tsit5", "abstol": 1e-8, "reltol": 1e-10}
    doses = [{"time": 0.0, "amount": 100.0}]

    # 1a. One-compartment IV
    results.append(run_benchmark(
        lambda: neopkpd.simulate(
            model_kind="OneCompIVBolus",
            params={"CL": 5.0, "V": 50.0},
            doses=doses,
            grid=grid,
            solver=solver
        ),
        "single_simulation", "simulation", "OneCompIV"
    ))

    # 1b. Two-compartment IV
    results.append(run_benchmark(
        lambda: neopkpd.simulate(
            model_kind="TwoCompIVBolus",
            params={"CL": 5.0, "Vc": 50.0, "Vp": 100.0, "Q": 10.0},
            doses=doses,
            grid=grid,
            solver=solver
        ),
        "single_simulation", "simulation", "TwoCompIV"
    ))

    # 1c. Two-compartment Oral
    results.append(run_benchmark(
        lambda: neopkpd.simulate(
            model_kind="TwoCompOral",
            params={"CL": 5.0, "Vc": 50.0, "Vp": 100.0, "Q": 10.0, "ka": 1.5},
            doses=doses,
            grid=grid,
            solver=solver
        ),
        "single_simulation", "simulation", "TwoCompOral"
    ))

    # 1d. Michaelis-Menten
    results.append(run_benchmark(
        lambda: neopkpd.simulate(
            model_kind="MichaelisMentenElimination",
            params={"Vmax": 50.0, "Km": 10.0, "V": 50.0},
            doses=doses,
            grid=grid,
            solver=solver
        ),
        "single_simulation", "simulation", "MichaelisMenten"
    ))

    return results


def benchmark_population_simulations():
    """Population simulation benchmarks."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Population Simulation Speed (Python)")
    print("=" * 60)

    results = []
    grid = {"tspan": (0.0, 24.0), "saveat": [i * 0.5 for i in range(49)]}
    solver = {"algorithm": "Tsit5", "abstol": 1e-8, "reltol": 1e-10}
    doses = [{"time": 0.0, "amount": 100.0}]

    for n_subj in [50, 100, 500, 1000]:
        results.append(run_benchmark(
            lambda n=n_subj: neopkpd.simulate_population(
                model_kind="TwoCompIVBolus",
                params={"CL": 5.0, "Vc": 50.0, "Vp": 100.0, "Q": 10.0},
                doses=doses,
                grid=grid,
                solver=solver,
                n_subjects=n,
                iiv={"CL": 0.3, "Vc": 0.2},
                seed=CONFIG["random_seed"]
            ),
            "population_simulation", "population", "TwoCompIV", n_subj
        ))

    return results


def benchmark_nca():
    """NCA computation benchmarks."""
    print("\n" + "=" * 60)
    print("BENCHMARK: NCA Computation Speed (Python)")
    print("=" * 60)

    results = []

    # Test data
    times = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0]
    conc = [0.0, 5.2, 8.1, 6.5, 4.2, 2.8, 1.9, 0.8, 0.1]

    from neopkpd.nca import nca_cmax, nca_tmax, nca_auc_0_t, nca_half_life

    # Individual NCA metrics
    results.append(run_benchmark(
        lambda: nca_cmax(conc),
        "nca_cmax", "nca", "Cmax"
    ))

    results.append(run_benchmark(
        lambda: nca_tmax(times, conc),
        "nca_tmax", "nca", "Tmax"
    ))

    results.append(run_benchmark(
        lambda: nca_auc_0_t(times, conc),
        "nca_auc", "nca", "AUC0t"
    ))

    results.append(run_benchmark(
        lambda: nca_half_life(times, conc),
        "nca_half_life", "nca", "HalfLife"
    ))

    return results


def benchmark_replay():
    """Artifact replay benchmarks."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Artifact Replay Speed (Python)")
    print("=" * 60)

    results = []

    # First, create an artifact by running a simulation
    grid = {"tspan": (0.0, 24.0), "saveat": [i * 0.5 for i in range(49)]}
    solver = {"algorithm": "Tsit5", "abstol": 1e-8, "reltol": 1e-10}
    doses = [{"time": 0.0, "amount": 100.0}]

    result = neopkpd.simulate(
        model_kind="TwoCompIVBolus",
        params={"CL": 5.0, "Vc": 50.0, "Vp": 100.0, "Q": 10.0},
        doses=doses,
        grid=grid,
        solver=solver,
        return_artifact=True
    )

    # Check if replay is available
    if hasattr(neopkpd, 'replay'):
        artifact_json = result.get('artifact', None)
        if artifact_json:
            results.append(run_benchmark(
                lambda: neopkpd.replay(artifact_json),
                "replay_artifact", "replay", "TwoCompIV"
            ))

    return results


# =============================================================================
# Main Runner
# =============================================================================

def run_all_benchmarks():
    """Run all benchmarks and save results."""
    all_results = []

    all_results.extend(benchmark_single_simulations())
    all_results.extend(benchmark_population_simulations())
    all_results.extend(benchmark_nca())
    all_results.extend(benchmark_replay())

    # Save results
    CONFIG["output_dir"].mkdir(parents=True, exist_ok=True)
    output_file = CONFIG["output_dir"] / "neopkpd_python_benchmarks.csv"

    with open(output_file, "w", newline="") as f:
        fieldnames = [
            "category", "name", "model", "n_subjects", "n_runs",
            "mean_ms", "std_ms", "median_ms", "min_ms", "max_ms",
            "ci_lower_ms", "ci_upper_ms", "timestamp",
            "python_version", "neopkpd_version"
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for result in all_results:
            writer.writerow(result.summarize())

    print("\n" + "=" * 70)
    print("BENCHMARK COMPLETE")
    print(f"Results saved to: {output_file}")
    print("=" * 70)

    return all_results


if __name__ == "__main__":
    run_all_benchmarks()
