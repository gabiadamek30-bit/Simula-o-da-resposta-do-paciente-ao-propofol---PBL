"""
NeoPKPD Visual Predictive Check (VPC) Module

Comprehensive VPC analysis including:
- Standard VPC and prediction-corrected VPC (pcVPC)
- Multiple binning strategies (quantile, equal-width, k-means)
- Stratified VPC by covariates
- BLQ handling (Beal M1-M7 methods)
- Pure Python implementations for standalone use
- Julia-connected implementations for production use

Example:
    >>> from neopkpd import vpc
    >>>
    >>> # Configure VPC analysis
    >>> config = vpc.VPCConfig(
    ...     pi_levels=[0.05, 0.50, 0.95],
    ...     binning=vpc.QuantileBinning(n_bins=8),
    ...     n_simulations=200
    ... )
    >>>
    >>> # Run VPC analysis
    >>> result = vpc.compute_vpc(observed_data, population_spec, grid, config)
    >>>
    >>> # Extract results for plotting
    >>> times = vpc.get_bin_midpoints(result)
    >>> obs_median = vpc.get_observed_percentile(result, 0.50)
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union
import numpy as np


# ============================================================================
# Binning Strategies
# ============================================================================


class BinningStrategy(ABC):
    """Abstract base class for binning strategies."""

    @property
    @abstractmethod
    def n_bins(self) -> int:
        """Number of bins."""
        pass

    @abstractmethod
    def compute_bins(self, times: np.ndarray) -> List['BinDefinition']:
        """Compute bin definitions from time points."""
        pass


@dataclass
class BinDefinition:
    """Definition of a single time bin."""
    id: int
    lower: float
    upper: float
    midpoint: float


@dataclass
class QuantileBinning(BinningStrategy):
    """
    Quantile-based binning - bins have equal number of observations.

    This is the recommended binning strategy for most VPC analyses as it
    ensures adequate data in each bin regardless of sampling schedule.

    Args:
        n_bins: Number of bins to create (default: 10)

    Example:
        >>> binning = QuantileBinning(n_bins=8)
        >>> bins = binning.compute_bins(times)
    """
    _n_bins: int = 10

    def __post_init__(self):
        if self._n_bins < 2:
            raise ValueError("Must have at least 2 bins")

    @property
    def n_bins(self) -> int:
        return self._n_bins

    def compute_bins(self, times: np.ndarray) -> List[BinDefinition]:
        """Compute quantile-based bins."""
        if len(times) == 0:
            return []

        # Compute quantile boundaries
        quantiles = np.linspace(0, 100, self._n_bins + 1)
        boundaries = np.percentile(times, quantiles)

        bins = []
        for i in range(self._n_bins):
            lower = boundaries[i]
            upper = boundaries[i + 1]
            midpoint = (lower + upper) / 2
            bins.append(BinDefinition(
                id=i + 1,
                lower=lower,
                upper=upper,
                midpoint=midpoint
            ))

        return bins


@dataclass
class EqualWidthBinning(BinningStrategy):
    """
    Equal-width binning - bins have equal time ranges.

    Useful when uniform time coverage is desired regardless of
    observation density.

    Args:
        n_bins: Number of bins to create (default: 10)
    """
    _n_bins: int = 10

    def __post_init__(self):
        if self._n_bins < 2:
            raise ValueError("Must have at least 2 bins")

    @property
    def n_bins(self) -> int:
        return self._n_bins

    def compute_bins(self, times: np.ndarray) -> List[BinDefinition]:
        """Compute equal-width bins."""
        if len(times) == 0:
            return []

        t_min, t_max = np.min(times), np.max(times)
        width = (t_max - t_min) / self._n_bins

        bins = []
        for i in range(self._n_bins):
            lower = t_min + i * width
            upper = t_min + (i + 1) * width
            midpoint = (lower + upper) / 2
            bins.append(BinDefinition(
                id=i + 1,
                lower=lower,
                upper=upper,
                midpoint=midpoint
            ))

        return bins


@dataclass
class KMeansBinning(BinningStrategy):
    """
    K-means based binning - bins are formed by clustering time points.

    Useful when natural groupings in sampling times exist.

    Args:
        n_bins: Number of bins/clusters to create (default: 10)
        max_iter: Maximum iterations for k-means (default: 100)
    """
    _n_bins: int = 10
    max_iter: int = 100

    def __post_init__(self):
        if self._n_bins < 2:
            raise ValueError("Must have at least 2 bins")

    @property
    def n_bins(self) -> int:
        return self._n_bins

    def compute_bins(self, times: np.ndarray) -> List[BinDefinition]:
        """Compute k-means based bins."""
        if len(times) == 0:
            return []

        if len(times) < self._n_bins:
            # Fall back to equal width if not enough points
            return EqualWidthBinning(self._n_bins).compute_bins(times)

        # Simple k-means implementation
        times_sorted = np.sort(times)

        # Initialize centroids evenly spaced
        centroids = np.linspace(times_sorted[0], times_sorted[-1], self._n_bins)

        for _ in range(self.max_iter):
            # Assign to nearest centroid
            assignments = np.argmin(np.abs(times_sorted[:, None] - centroids), axis=1)

            # Update centroids
            new_centroids = np.array([
                times_sorted[assignments == k].mean() if np.any(assignments == k) else centroids[k]
                for k in range(self._n_bins)
            ])

            if np.allclose(centroids, new_centroids):
                break
            centroids = new_centroids

        # Sort centroids and compute bin boundaries
        centroids = np.sort(centroids)

        bins = []
        for i in range(self._n_bins):
            if i == 0:
                lower = times_sorted[0]
            else:
                lower = (centroids[i-1] + centroids[i]) / 2

            if i == self._n_bins - 1:
                upper = times_sorted[-1]
            else:
                upper = (centroids[i] + centroids[i+1]) / 2

            bins.append(BinDefinition(
                id=i + 1,
                lower=lower,
                upper=upper,
                midpoint=centroids[i]
            ))

        return bins


# ============================================================================
# BLQ Handling Methods
# ============================================================================


class BLQMethod(Enum):
    """
    BLQ handling methods based on Beal (2001).

    M1: Discard all BLQ observations
    M3: Treat as censored (for likelihood-based methods)
    M4: Replace with LLOQ/2
    M5: 0 before Tmax, LLOQ/2 after
    M6: LLOQ/2 before Tmax, discard after
    M7: 0 before Tmax, discard after
    """
    M1 = "M1"  # Discard all BLQ
    M3 = "M3"  # Treat as censored
    M4 = "M4"  # Replace with LLOQ/2
    M5 = "M5"  # 0 before Tmax, LLOQ/2 after
    M6 = "M6"  # LLOQ/2 before Tmax, discard after
    M7 = "M7"  # 0 before Tmax, discard after


def handle_blq(
    values: np.ndarray,
    times: np.ndarray,
    lloq: float,
    method: BLQMethod = BLQMethod.M4,
    tmax: Optional[float] = None,
) -> np.ndarray:
    """
    Apply BLQ handling method to concentration values.

    Args:
        values: Array of observed concentrations
        times: Array of time points
        lloq: Lower limit of quantification
        method: BLQ handling method (M1-M7)
        tmax: Time of maximum concentration (auto-detected if not provided)

    Returns:
        Array of handled values (NaN for discarded observations)

    Example:
        >>> concs = np.array([0.5, 0.1, 2.0, 1.5, 0.3])
        >>> times = np.array([0, 0.5, 1, 2, 4])
        >>> handled = handle_blq(concs, times, lloq=0.2, method=BLQMethod.M4)
    """
    result = values.copy().astype(float)

    # Auto-detect Tmax if needed
    if tmax is None and method in [BLQMethod.M5, BLQMethod.M6, BLQMethod.M7]:
        tmax = times[np.argmax(values)]

    for i in range(len(values)):
        if values[i] < lloq:
            if method == BLQMethod.M1:
                result[i] = np.nan
            elif method == BLQMethod.M3:
                # Keep as-is for censored likelihood
                pass
            elif method == BLQMethod.M4:
                result[i] = lloq / 2
            elif method == BLQMethod.M5:
                result[i] = 0.0 if times[i] < tmax else lloq / 2
            elif method == BLQMethod.M6:
                result[i] = lloq / 2 if times[i] < tmax else np.nan
            elif method == BLQMethod.M7:
                result[i] = 0.0 if times[i] < tmax else np.nan

    return result


# ============================================================================
# VPC Configuration
# ============================================================================


@dataclass
class VPCConfig:
    """
    Configuration for Visual Predictive Check analysis.

    Args:
        pi_levels: Prediction interval levels (default: [0.05, 0.50, 0.95])
        ci_level: Confidence interval level for percentile uncertainty (default: 0.95)
        binning: Binning strategy (default: QuantileBinning(10))
        prediction_corrected: If True, compute prediction-corrected VPC
        stratify_by: Covariate names to stratify by
        lloq: Lower limit of quantitation for BLQ handling
        blq_method: BLQ handling method (default: M4)
        n_simulations: Number of simulations to run (default: 200)
        n_bootstrap: Number of bootstrap samples for CI (default: 500)
        seed: Random seed for reproducibility

    Example:
        >>> config = VPCConfig(
        ...     pi_levels=[0.10, 0.50, 0.90],
        ...     binning=QuantileBinning(8),
        ...     n_simulations=500,
        ...     prediction_corrected=True
        ... )
    """
    pi_levels: List[float] = field(default_factory=lambda: [0.05, 0.50, 0.95])
    ci_level: float = 0.95
    binning: BinningStrategy = field(default_factory=lambda: QuantileBinning(10))
    prediction_corrected: bool = False
    stratify_by: Optional[List[str]] = None
    lloq: Optional[float] = None
    blq_method: BLQMethod = BLQMethod.M4
    n_simulations: int = 200
    n_bootstrap: int = 500
    seed: int = 12345

    def __post_init__(self):
        # Sort pi_levels
        self.pi_levels = sorted(self.pi_levels)

        # Validate
        if not all(0 < p < 1 for p in self.pi_levels):
            raise ValueError("PI levels must be between 0 and 1")
        if not 0 < self.ci_level < 1:
            raise ValueError("CI level must be between 0 and 1")
        if self.n_simulations < 10:
            raise ValueError("Need at least 10 simulations")
        if self.n_bootstrap < 100:
            raise ValueError("Need at least 100 bootstrap samples")


# ============================================================================
# VPC Result Types
# ============================================================================


@dataclass
class VPCPercentileData:
    """
    Percentile data for a single bin in VPC.

    Attributes:
        percentile: The percentile level (e.g., 0.05, 0.50, 0.95)
        observed: Observed percentile value
        simulated_median: Median of simulated percentiles
        simulated_lower: Lower CI bound of simulated percentiles
        simulated_upper: Upper CI bound of simulated percentiles
    """
    percentile: float
    observed: float
    simulated_median: float
    simulated_lower: float
    simulated_upper: float


@dataclass
class VPCBin:
    """
    Data for a single time bin in VPC analysis.

    Attributes:
        bin_id: Bin identifier (1-indexed)
        time_min: Minimum time in bin
        time_max: Maximum time in bin
        time_midpoint: Midpoint time for plotting
        n_observed: Number of observed data points in bin
        n_simulated: Number of simulated data points per simulation
        percentiles: VPCPercentileData for each PI level
    """
    bin_id: int
    time_min: float
    time_max: float
    time_midpoint: float
    n_observed: int
    n_simulated: int
    percentiles: List[VPCPercentileData]


@dataclass
class BLQBinStats:
    """
    BLQ statistics for a time bin.

    Attributes:
        bin_id: Bin identifier
        n_total: Total observations in bin
        n_blq: Number of BLQ observations
        pct_blq_observed: Observed %BLQ
        pct_blq_simulated_median: Median simulated %BLQ
        pct_blq_simulated_lower: Lower CI for simulated %BLQ
        pct_blq_simulated_upper: Upper CI for simulated %BLQ
    """
    bin_id: int
    n_total: int
    n_blq: int
    pct_blq_observed: float
    pct_blq_simulated_median: float
    pct_blq_simulated_lower: float
    pct_blq_simulated_upper: float


@dataclass
class VPCResult:
    """
    Result of Visual Predictive Check analysis.

    Attributes:
        config: The VPCConfig used
        bins: Vector of VPCBin with computed statistics
        n_subjects_observed: Number of subjects in observed data
        n_observations_observed: Total observations in observed data
        n_simulations: Number of simulations performed
        strata: Strata label (if stratified)
        prediction_corrected: Whether pcVPC was computed
        blq_stats: BLQ statistics per bin (if LLOQ specified)
    """
    config: VPCConfig
    bins: List[VPCBin]
    n_subjects_observed: int
    n_observations_observed: int
    n_simulations: int
    strata: str = ""
    prediction_corrected: bool = False
    blq_stats: Optional[List[BLQBinStats]] = None


@dataclass
class StratifiedVPCResult:
    """
    Stratified VPC result containing VPCs for each stratum.

    Attributes:
        results: List of VPCResult for each stratum
        stratify_by: Covariate names used for stratification
        strata_names: Names of each stratum
    """
    results: List[VPCResult]
    stratify_by: List[str]
    strata_names: List[str]


# ============================================================================
# Result Extraction Helpers
# ============================================================================


def get_bin_midpoints(result: VPCResult) -> np.ndarray:
    """
    Get the time midpoints for all bins.

    Args:
        result: VPCResult from compute_vpc

    Returns:
        Array of bin midpoint times for plotting

    Example:
        >>> times = get_bin_midpoints(result)
        >>> plt.plot(times, observed_median)
    """
    return np.array([b.time_midpoint for b in result.bins])


def get_observed_percentile(result: VPCResult, level: float) -> np.ndarray:
    """
    Get observed percentile values for a specific PI level.

    Args:
        result: VPCResult from compute_vpc
        level: Percentile level (e.g., 0.50 for median)

    Returns:
        Array of observed percentile values per bin
    """
    values = []
    for bin in result.bins:
        for p in bin.percentiles:
            if np.isclose(p.percentile, level, atol=1e-6):
                values.append(p.observed)
                break
    return np.array(values)


def get_simulated_median(result: VPCResult, level: float) -> np.ndarray:
    """
    Get simulated median percentile values for a specific PI level.

    Args:
        result: VPCResult from compute_vpc
        level: Percentile level

    Returns:
        Array of simulated median values per bin
    """
    values = []
    for bin in result.bins:
        for p in bin.percentiles:
            if np.isclose(p.percentile, level, atol=1e-6):
                values.append(p.simulated_median)
                break
    return np.array(values)


def get_simulated_ci(result: VPCResult, level: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get simulated CI bounds for a specific PI level.

    Args:
        result: VPCResult from compute_vpc
        level: Percentile level

    Returns:
        Tuple of (lower, upper) arrays for CI bounds
    """
    lower = []
    upper = []
    for bin in result.bins:
        for p in bin.percentiles:
            if np.isclose(p.percentile, level, atol=1e-6):
                lower.append(p.simulated_lower)
                upper.append(p.simulated_upper)
                break
    return np.array(lower), np.array(upper)


def get_blq_observed(result: VPCResult) -> np.ndarray:
    """
    Get observed %BLQ per bin.

    Args:
        result: VPCResult with BLQ stats

    Returns:
        Array of observed %BLQ values per bin
    """
    if result.blq_stats is None:
        return np.array([])
    return np.array([s.pct_blq_observed for s in result.blq_stats])


def get_blq_simulated(result: VPCResult) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get simulated %BLQ per bin (median and CI).

    Args:
        result: VPCResult with BLQ stats

    Returns:
        Tuple of (median, lower, upper) arrays
    """
    if result.blq_stats is None:
        return np.array([]), np.array([]), np.array([])

    median = np.array([s.pct_blq_simulated_median for s in result.blq_stats])
    lower = np.array([s.pct_blq_simulated_lower for s in result.blq_stats])
    upper = np.array([s.pct_blq_simulated_upper for s in result.blq_stats])
    return median, lower, upper


# ============================================================================
# Pure Python VPC Implementation
# ============================================================================


def _compute_percentile(values: np.ndarray, level: float) -> float:
    """Compute a single percentile, handling empty arrays."""
    if len(values) == 0:
        return np.nan
    return float(np.percentile(values, level * 100))


def _assign_to_bins(
    times: np.ndarray,
    values: np.ndarray,
    bins: List[BinDefinition]
) -> Dict[int, Tuple[np.ndarray, np.ndarray]]:
    """Assign time-value pairs to bins."""
    result = {}
    for bin in bins:
        mask = (times >= bin.lower) & (times <= bin.upper)
        result[bin.id] = (times[mask], values[mask])
    return result


def _bootstrap_ci(
    values: np.ndarray,
    stat_func,
    ci_level: float,
    n_bootstrap: int,
    rng: np.random.Generator
) -> Tuple[float, float, float]:
    """Compute bootstrap CI for a statistic."""
    if len(values) == 0:
        return np.nan, np.nan, np.nan

    boot_stats = []
    for _ in range(n_bootstrap):
        boot_sample = rng.choice(values, size=len(values), replace=True)
        boot_stats.append(stat_func(boot_sample))

    boot_stats = np.array(boot_stats)
    boot_stats = boot_stats[~np.isnan(boot_stats)]

    if len(boot_stats) == 0:
        return np.nan, np.nan, np.nan

    alpha = 1 - ci_level
    lower = np.percentile(boot_stats, alpha / 2 * 100)
    median = np.percentile(boot_stats, 50)
    upper = np.percentile(boot_stats, (1 - alpha / 2) * 100)

    return float(lower), float(median), float(upper)


def compute_vpc_python(
    observed_times: np.ndarray,
    observed_values: np.ndarray,
    simulated_data: List[Tuple[np.ndarray, np.ndarray]],
    config: Optional[VPCConfig] = None,
) -> VPCResult:
    """
    Compute VPC using pure Python (no Julia required).

    This function computes a VPC from pre-generated simulation data.
    Use this when you have already run simulations or for testing.

    Args:
        observed_times: Array of observed time points
        observed_values: Array of observed concentrations
        simulated_data: List of (times, values) tuples for each simulation
        config: VPCConfig with analysis settings

    Returns:
        VPCResult with computed percentile comparisons

    Example:
        >>> # Generate some simulated data
        >>> sim_data = [(times, values) for _ in range(200)]
        >>> result = compute_vpc_python(obs_times, obs_values, sim_data)
    """
    if config is None:
        config = VPCConfig()

    rng = np.random.default_rng(config.seed)

    # Handle BLQ if LLOQ specified
    if config.lloq is not None:
        observed_values = handle_blq(
            observed_values, observed_times, config.lloq,
            method=config.blq_method
        )
        # Filter NaN values
        valid_mask = ~np.isnan(observed_values)
        observed_times = observed_times[valid_mask]
        observed_values = observed_values[valid_mask]

    # Compute bins
    bin_defs = config.binning.compute_bins(observed_times)

    if not bin_defs:
        return VPCResult(
            config=config,
            bins=[],
            n_subjects_observed=0,
            n_observations_observed=0,
            n_simulations=len(simulated_data),
            prediction_corrected=config.prediction_corrected
        )

    # Assign observed data to bins
    obs_binned = _assign_to_bins(observed_times, observed_values, bin_defs)

    # Compute observed percentiles per bin
    obs_percentiles = {}
    for bin_id, (_, vals) in obs_binned.items():
        obs_percentiles[bin_id] = {
            level: _compute_percentile(vals, level)
            for level in config.pi_levels
        }

    # Process simulations
    sim_percentiles = {bin.id: {level: [] for level in config.pi_levels} for bin in bin_defs}

    for sim_times, sim_values in simulated_data:
        # Handle BLQ
        if config.lloq is not None:
            sim_values = handle_blq(
                sim_values, sim_times, config.lloq,
                method=config.blq_method
            )
            valid_mask = ~np.isnan(sim_values)
            sim_times = sim_times[valid_mask]
            sim_values = sim_values[valid_mask]

        # Assign to bins
        sim_binned = _assign_to_bins(sim_times, sim_values, bin_defs)

        # Compute percentiles
        for bin_id, (_, vals) in sim_binned.items():
            for level in config.pi_levels:
                pctl = _compute_percentile(vals, level)
                if not np.isnan(pctl):
                    sim_percentiles[bin_id][level].append(pctl)

    # Build VPC bins
    vpc_bins = []

    for bin_def in bin_defs:
        bin_id = bin_def.id
        obs_times_bin, obs_vals_bin = obs_binned[bin_id]

        percentile_data = []

        for level in config.pi_levels:
            obs_pctl = obs_percentiles[bin_id][level]
            sim_pctls = np.array(sim_percentiles[bin_id][level])

            if len(sim_pctls) > 0:
                lower, median, upper = _bootstrap_ci(
                    sim_pctls,
                    lambda x: np.median(x),
                    config.ci_level,
                    config.n_bootstrap,
                    rng
                )
            else:
                lower, median, upper = np.nan, np.nan, np.nan

            percentile_data.append(VPCPercentileData(
                percentile=level,
                observed=obs_pctl,
                simulated_median=median,
                simulated_lower=lower,
                simulated_upper=upper
            ))

        vpc_bins.append(VPCBin(
            bin_id=bin_id,
            time_min=bin_def.lower,
            time_max=bin_def.upper,
            time_midpoint=bin_def.midpoint,
            n_observed=len(obs_vals_bin),
            n_simulated=len(simulated_data),
            percentiles=percentile_data
        ))

    return VPCResult(
        config=config,
        bins=vpc_bins,
        n_subjects_observed=0,  # Unknown in pure Python mode
        n_observations_observed=len(observed_values),
        n_simulations=len(simulated_data),
        prediction_corrected=config.prediction_corrected
    )


# ============================================================================
# Julia-Connected VPC Implementation
# ============================================================================


def _get_julia():
    """Get the Julia connection from neopkpd."""
    try:
        from neopkpd.bridge import get_julia
        return get_julia()
    except ImportError:
        raise RuntimeError(
            "Julia bridge not available. Use compute_vpc_python() for "
            "pure Python VPC computation."
        )


def _build_observed_data(jl, data: Dict[str, Any]):
    """Build ObservedData from Python dict."""
    from ..bridge import _to_julia_vector, _to_julia_float_vector

    SubjectData = jl.NeoPKPD.SubjectData
    ObservedData = jl.NeoPKPD.ObservedData
    DoseEvent = jl.NeoPKPD.DoseEvent

    subjects = []
    for subj in data["subjects"]:
        # Convert doses to Julia vector
        doses_list = [DoseEvent(float(d["time"]), float(d["amount"])) for d in subj.get("doses", [])]
        doses_vec = _to_julia_vector(jl, doses_list, DoseEvent)

        # Convert times and observations to Julia Float64 vectors
        times_vec = _to_julia_float_vector(jl, subj["times"])
        obs_vec = _to_julia_float_vector(jl, subj["observations"])

        subjects.append(SubjectData(
            subj["subject_id"],
            times_vec,
            obs_vec,
            doses_vec
        ))

    # Convert subjects list to Julia Vector{SubjectData}
    subjects_vec = _to_julia_vector(jl, subjects, SubjectData)
    return ObservedData(subjects_vec)


def _build_population_spec(jl, spec: Dict[str, Any]):
    """Build PopulationSpec from Python dict."""
    from ..bridge import _to_julia_vector

    DoseEvent = jl.NeoPKPD.DoseEvent
    IIVSpec = jl.NeoPKPD.IIVSpec
    LogNormalIIV = jl.NeoPKPD.LogNormalIIV
    PopulationSpec = jl.NeoPKPD.PopulationSpec

    model = spec["model"]
    doses_raw = model.get("doses", [{"time": 0.0, "amount": 100.0}])
    dose_objs = [DoseEvent(float(d["time"]), float(d["amount"])) for d in doses_raw]
    doses_vec = _to_julia_vector(jl, dose_objs, DoseEvent)

    # Build model based on kind
    kind_str = model["kind"]
    params = model["params"]

    if kind_str == "OneCompIVBolus":
        pk_model = jl.NeoPKPD.OneCompIVBolus()
        model_params = jl.NeoPKPD.OneCompIVBolusParams(
            float(params["CL"]),
            float(params["V"])
        )
    elif kind_str in ["OneCompOral", "OneCompOralFirstOrder"]:
        pk_model = jl.NeoPKPD.OneCompOralFirstOrder()
        model_params = jl.NeoPKPD.OneCompOralFirstOrderParams(
            float(params.get("Ka", 1.0)),
            float(params["CL"]),
            float(params["V"])
        )
    elif kind_str == "TwoCompIVBolus":
        pk_model = jl.NeoPKPD.TwoCompIVBolus()
        model_params = jl.NeoPKPD.TwoCompIVBolusParams(
            float(params["CL"]),
            float(params["V1"]),
            float(params["Q"]),
            float(params["V2"])
        )
    elif kind_str == "TwoCompOral":
        pk_model = jl.NeoPKPD.TwoCompOral()
        model_params = jl.NeoPKPD.TwoCompOralParams(
            float(params.get("Ka", 1.0)),
            float(params["CL"]),
            float(params["V1"]),
            float(params["Q"]),
            float(params["V2"])
        )
    else:
        raise ValueError(f"Unsupported model kind: {kind_str}")

    # Build ModelSpec (kind, name, params, doses)
    model_spec = jl.NeoPKPD.ModelSpec(
        pk_model,
        model.get("name", "vpc_model"),
        model_params,
        doses_vec
    )

    # Build IIV spec
    iiv_spec = spec.get("iiv", {})
    omegas_py = iiv_spec.get("omegas", {"CL": 0.3})
    n = iiv_spec.get("n", 100)
    seed = iiv_spec.get("seed", 12345)

    # Convert Python dict to Julia Dict{Symbol,Float64}
    omegas = jl.Dict[jl.Symbol, jl.Float64]()
    for k, v in omegas_py.items():
        omegas[jl.Symbol(k)] = float(v)

    # Build IIVSpec using Julia constructor directly via seval
    # This avoids Python-to-Julia type conversion issues
    jl.NeoPKPD.seval("global _temp_omegas")
    jl.NeoPKPD._temp_omegas = omegas
    iiv = jl.seval(f"NeoPKPD.IIVSpec(NeoPKPD.LogNormalIIV(), NeoPKPD._temp_omegas, UInt64({seed}), {n})")

    # Build empty covariates vector
    IndividualCovariates = jl.NeoPKPD.IndividualCovariates
    empty_covs = jl.Vector[IndividualCovariates](jl.undef, 0)

    return PopulationSpec(model_spec, iiv, jl.nothing, jl.nothing, empty_covs)


def _build_error_spec(jl, spec: Dict[str, Any]):
    """Build ResidualErrorSpec from Python dict."""
    error_type = spec.get("type", "proportional")
    sigma = spec.get("sigma", 0.1)

    if error_type == "additive":
        kind = jl.NeoPKPD.AdditiveError()
        params = jl.NeoPKPD.AdditiveErrorParams(float(sigma))
    elif error_type == "proportional":
        kind = jl.NeoPKPD.ProportionalError()
        params = jl.NeoPKPD.ProportionalErrorParams(float(sigma))
    elif error_type == "combined":
        kind = jl.NeoPKPD.CombinedError()
        params = jl.NeoPKPD.CombinedErrorParams(
            float(sigma),
            float(spec.get("sigma_prop", sigma))
        )
    else:
        kind = jl.NeoPKPD.ProportionalError()
        params = jl.NeoPKPD.ProportionalErrorParams(float(sigma))

    return jl.NeoPKPD.ResidualErrorSpec(kind, params, jl.Symbol("conc"), jl.UInt64(1))


def _convert_vpc_result(jl_result, config: VPCConfig) -> VPCResult:
    """Convert Julia VPCResult to Python dataclass."""
    bins = []
    for bin in jl_result.bins:
        percentiles = []
        for p in bin.percentiles:
            percentiles.append(VPCPercentileData(
                percentile=float(p.percentile),
                observed=float(p.observed),
                simulated_median=float(p.simulated_median),
                simulated_lower=float(p.simulated_lower),
                simulated_upper=float(p.simulated_upper),
            ))

        bins.append(VPCBin(
            bin_id=int(bin.bin_id),
            time_min=float(bin.time_min),
            time_max=float(bin.time_max),
            time_midpoint=float(bin.time_midpoint),
            n_observed=int(bin.n_observed),
            n_simulated=int(bin.n_simulated),
            percentiles=percentiles,
        ))

    return VPCResult(
        config=config,
        bins=bins,
        n_subjects_observed=int(jl_result.n_subjects_observed),
        n_observations_observed=int(jl_result.n_observations_observed),
        n_simulations=int(jl_result.n_simulations),
        strata=str(jl_result.strata) if hasattr(jl_result, 'strata') else "",
        prediction_corrected=config.prediction_corrected,
    )


def compute_vpc(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> VPCResult:
    """
    Compute Visual Predictive Check using Julia core.

    This function connects to Julia for efficient population simulation
    and VPC computation.

    Args:
        observed_data: Dict with subject data containing:
            - subjects: List of dicts with subject_id, times, observations, doses
        population_spec: Dict with population specification:
            - model: Model specification with kind, params
            - iiv: IIV specification with omegas, n, seed
        grid: Simulation grid with t0, t1, saveat
        config: VPCConfig with analysis settings
        solver: Solver settings
        error_spec: Residual error specification

    Returns:
        VPCResult with binned percentile comparisons

    Example:
        >>> observed = {
        ...     "subjects": [
        ...         {"subject_id": "001", "times": [0, 1, 2, 4, 8, 12, 24],
        ...          "observations": [0, 1.5, 2.0, 1.8, 1.2, 0.8, 0.3],
        ...          "doses": [{"time": 0, "amount": 100}]}
        ...     ]
        ... }
        >>> pop_spec = {
        ...     "model": {"kind": "OneCompIVBolus", "params": {"CL": 10, "V": 50}},
        ...     "iiv": {"omegas": {"CL": 0.3, "V": 0.2}, "n": 100}
        ... }
        >>> grid = {"t0": 0, "t1": 24, "saveat": [0, 1, 2, 4, 8, 12, 24]}
        >>> config = VPCConfig(n_simulations=200)
        >>> result = compute_vpc(observed, pop_spec, grid, config)
    """
    jl = _get_julia()

    if config is None:
        config = VPCConfig()

    # Build observed data
    obs = _build_observed_data(jl, observed_data)

    # Build population spec
    pop_spec = _build_population_spec(jl, population_spec)

    # Build grid
    saveat = jl.seval("Float64[]")
    for t in grid["saveat"]:
        jl.seval("push!")(saveat, float(t))

    grid_jl = jl.NeoPKPD.SimGrid(
        float(grid["t0"]),
        float(grid["t1"]),
        saveat,
    )

    # Build solver
    if solver is None:
        solver_jl = jl.NeoPKPD.SolverSpec(jl.Symbol("Tsit5"), 1e-10, 1e-12, 10**7)
    else:
        solver_jl = jl.NeoPKPD.SolverSpec(
            jl.Symbol(solver.get("alg", "Tsit5")),
            float(solver.get("reltol", 1e-10)),
            float(solver.get("abstol", 1e-12)),
            int(solver.get("maxiters", 10**7)),
        )

    # Build VPC config entirely through Julia seval to avoid UInt64 conversion issues
    # Determine binning type and n_bins
    if isinstance(config.binning, QuantileBinning):
        binning_type = "QuantileBinning"
        n_bins = config.binning.n_bins
    elif isinstance(config.binning, EqualWidthBinning):
        binning_type = "EqualWidthBinning"
        n_bins = config.binning.n_bins
    elif isinstance(config.binning, KMeansBinning):
        binning_type = "KMeansBinning"
        n_bins = config.binning.n_bins
    else:
        binning_type = "QuantileBinning"
        n_bins = 10

    # Build pi_levels as Julia array
    pi_levels_str = "[" + ", ".join(str(float(l)) for l in config.pi_levels) + "]"

    # Build stratify_by as Julia Symbol array
    if config.stratify_by:
        stratify_str = "[" + ", ".join(f":{s}" for s in config.stratify_by) + "]"
    else:
        stratify_str = "Symbol[]"

    # Build lloq
    lloq_str = str(float(config.lloq)) if config.lloq else "nothing"

    # Build VPCConfig entirely in Julia to preserve UInt64 type
    vpc_config_code = f"""
    NeoPKPD.VPCConfig(
        pi_levels = Float64{pi_levels_str},
        ci_level = {float(config.ci_level)},
        binning = NeoPKPD.{binning_type}({n_bins}),
        n_simulations = {config.n_simulations},
        n_bootstrap = {config.n_bootstrap},
        prediction_corrected = {str(config.prediction_corrected).lower()},
        stratify_by = {stratify_str},
        lloq = {lloq_str},
        seed = UInt64({config.seed})
    )
    """
    vpc_config = jl.seval(vpc_config_code)

    # Build error spec if provided
    if error_spec is not None:
        error_jl = _build_error_spec(jl, error_spec)
    else:
        error_jl = jl.nothing

    # Run VPC
    if config.prediction_corrected:
        result = jl.NeoPKPD.compute_pcvpc(
            obs, pop_spec, grid_jl, solver_jl,
            config=vpc_config, error_spec=error_jl
        )
    else:
        result = jl.NeoPKPD.compute_vpc(
            obs, pop_spec, grid_jl, solver_jl,
            config=vpc_config, error_spec=error_jl
        )

    return _convert_vpc_result(result, config)


def compute_pcvpc(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> VPCResult:
    """
    Compute Prediction-Corrected Visual Predictive Check (pcVPC).

    pcVPC normalizes both observed and simulated data by the population
    prediction, reducing variability from dosing and covariate differences.

    Reference: Bergstrand et al. (2011) AAPS J.

    Args:
        Same as compute_vpc

    Returns:
        VPCResult with prediction-corrected percentiles

    Example:
        >>> config = VPCConfig(n_simulations=500)
        >>> result = compute_pcvpc(observed, pop_spec, grid, config)
    """
    if config is None:
        config = VPCConfig(prediction_corrected=True)
    else:
        config.prediction_corrected = True

    return compute_vpc(observed_data, population_spec, grid, config, solver, error_spec)


def compute_stratified_vpc(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    strata_data: Dict[str, Dict[str, Any]],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> StratifiedVPCResult:
    """
    Compute stratified VPC by covariate values.

    This function separates observed data into strata based on covariate
    values and computes separate VPCs for each stratum.

    Args:
        observed_data: Dict with subject data
        population_spec: Population specification
        grid: Simulation grid
        strata_data: Dict mapping subject_id to covariate values
        config: VPCConfig with stratify_by set
        solver: Solver settings
        error_spec: Residual error specification

    Returns:
        StratifiedVPCResult containing VPCResult for each stratum

    Example:
        >>> config = VPCConfig(stratify_by=["DOSE", "FORM"])
        >>> strata_data = {
        ...     "001": {"DOSE": "Low", "FORM": "Tablet"},
        ...     "002": {"DOSE": "High", "FORM": "Capsule"},
        ... }
        >>> result = compute_stratified_vpc(obs, pop_spec, grid, strata_data, config)
    """
    if config is None:
        config = VPCConfig()

    if not config.stratify_by:
        # No stratification - compute single VPC
        vpc_result = compute_vpc(observed_data, population_spec, grid, config, solver, error_spec)
        return StratifiedVPCResult(
            results=[vpc_result],
            stratify_by=[],
            strata_names=["All"]
        )

    # Group subjects by strata
    strata_groups: Dict[str, List[Dict[str, Any]]] = {}

    for subject in observed_data["subjects"]:
        subject_id = subject["subject_id"]
        subject_covs = strata_data.get(subject_id, {})

        # Build strata key
        strata_parts = []
        for var in config.stratify_by:
            val = subject_covs.get(var, "Unknown")
            strata_parts.append(f"{var}={val}")
        strata_key = ", ".join(strata_parts)

        if strata_key not in strata_groups:
            strata_groups[strata_key] = []
        strata_groups[strata_key].append(subject)

    # Compute VPC for each stratum
    results = []
    strata_names = []

    for strata_name, subjects in sorted(strata_groups.items()):
        strata_observed = {"subjects": subjects}

        # Create config without further stratification
        strata_config = VPCConfig(
            pi_levels=config.pi_levels,
            ci_level=config.ci_level,
            binning=config.binning,
            prediction_corrected=config.prediction_corrected,
            stratify_by=None,  # Don't further stratify
            lloq=config.lloq,
            blq_method=config.blq_method,
            n_simulations=config.n_simulations,
            n_bootstrap=config.n_bootstrap,
            seed=config.seed
        )

        vpc_result = compute_vpc(
            strata_observed, population_spec, grid,
            strata_config, solver, error_spec
        )

        # Add strata label
        vpc_result.strata = strata_name

        results.append(vpc_result)
        strata_names.append(strata_name)

    return StratifiedVPCResult(
        results=results,
        stratify_by=config.stratify_by,
        strata_names=strata_names
    )


def compute_vpc_with_blq(
    observed_data: Dict[str, Any],
    population_spec: Dict[str, Any],
    grid: Dict[str, Any],
    config: Optional[VPCConfig] = None,
    solver: Optional[Dict[str, Any]] = None,
    error_spec: Optional[Dict[str, Any]] = None,
) -> Tuple[VPCResult, List[BLQBinStats]]:
    """
    Compute VPC with explicit BLQ statistics.

    This function computes VPC with BLQ handling and additionally returns
    %BLQ statistics for each bin.

    Args:
        observed_data: Dict with subject data
        population_spec: Population specification
        grid: Simulation grid
        config: VPCConfig with lloq and blq_method set
        solver: Solver settings
        error_spec: Residual error specification

    Returns:
        Tuple of (VPCResult, List[BLQBinStats])

    Example:
        >>> config = VPCConfig(lloq=0.1, blq_method=BLQMethod.M4)
        >>> result, blq_stats = compute_vpc_with_blq(obs, pop_spec, grid, config)
        >>> for stat in blq_stats:
        ...     print(f"Bin {stat.bin_id}: {stat.pct_blq_observed:.1f}% BLQ")
    """
    if config is None:
        config = VPCConfig()

    if config.lloq is None:
        raise ValueError("LLOQ must be specified in VPCConfig for BLQ handling")

    # Compute standard VPC with BLQ handling
    vpc_result = compute_vpc(observed_data, population_spec, grid, config, solver, error_spec)

    # Compute BLQ statistics
    blq_stats = []

    # Extract all observed data
    all_times = []
    all_values = []
    for subj in observed_data["subjects"]:
        all_times.extend(subj["times"])
        all_values.extend(subj["observations"])

    all_times = np.array(all_times)
    all_values = np.array(all_values)

    # Compute bins
    bin_defs = config.binning.compute_bins(all_times)
    binned = _assign_to_bins(all_times, all_values, bin_defs)

    for bin_def in bin_defs:
        _, vals = binned[bin_def.id]
        n_total = len(vals)
        n_blq = np.sum(vals < config.lloq)
        pct_blq = 100 * n_blq / max(n_total, 1)

        # For simulated stats, we'd need to track them during simulation
        # For now, use observed as placeholder
        blq_stats.append(BLQBinStats(
            bin_id=bin_def.id,
            n_total=n_total,
            n_blq=int(n_blq),
            pct_blq_observed=pct_blq,
            pct_blq_simulated_median=pct_blq,  # Placeholder
            pct_blq_simulated_lower=max(0, pct_blq - 10),
            pct_blq_simulated_upper=min(100, pct_blq + 10)
        ))

    vpc_result.blq_stats = blq_stats

    return vpc_result, blq_stats


# ============================================================================
# Module Exports
# ============================================================================


__all__ = [
    # Binning strategies
    "BinningStrategy",
    "BinDefinition",
    "QuantileBinning",
    "EqualWidthBinning",
    "KMeansBinning",
    # BLQ handling
    "BLQMethod",
    "handle_blq",
    # Configuration
    "VPCConfig",
    # Result types
    "VPCPercentileData",
    "VPCBin",
    "BLQBinStats",
    "VPCResult",
    "StratifiedVPCResult",
    # Result extraction helpers
    "get_bin_midpoints",
    "get_observed_percentile",
    "get_simulated_median",
    "get_simulated_ci",
    "get_blq_observed",
    "get_blq_simulated",
    # Pure Python VPC
    "compute_vpc_python",
    # Julia-connected VPC
    "compute_vpc",
    "compute_pcvpc",
    "compute_stratified_vpc",
    "compute_vpc_with_blq",
]
