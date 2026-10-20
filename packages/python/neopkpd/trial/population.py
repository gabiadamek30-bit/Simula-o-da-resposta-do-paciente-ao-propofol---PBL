"""
NeoPKPD Virtual Population Generation

Virtual population generation for clinical trial simulation including:
- Demographic specifications
- Disease state modeling
- Covariate correlations
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
import random
import math


@dataclass
class DemographicSpec:
    """
    Demographic distribution specification.

    Attributes:
        age_mean: Mean age (years)
        age_sd: Age standard deviation
        age_range: Age range limits (min, max)
        weight_mean: Mean weight (kg)
        weight_sd: Weight standard deviation
        weight_range: Weight range limits
        female_proportion: Proportion of females (0-1)
        race_distribution: Race/ethnicity distribution

    Example:
        >>> spec = DemographicSpec(age_mean=45.0, age_sd=12.0)
    """
    age_mean: float = 45.0
    age_sd: float = 15.0
    age_range: Tuple[float, float] = (18.0, 75.0)
    weight_mean: float = 75.0
    weight_sd: float = 15.0
    weight_range: Tuple[float, float] = (40.0, 150.0)
    female_proportion: float = 0.5
    race_distribution: Dict[str, float] = field(default_factory=lambda: {
        "caucasian": 0.70,
        "asian": 0.15,
        "black": 0.10,
        "other": 0.05,
    })


@dataclass
class DiseaseSpec:
    """
    Disease state specification.

    Attributes:
        name: Disease name
        severity_distribution: Distribution of severity levels
        baseline_biomarker_mean: Mean baseline biomarker value
        baseline_biomarker_sd: Baseline biomarker standard deviation

    Example:
        >>> spec = DiseaseSpec(
        ...     name="diabetes",
        ...     severity_distribution={"mild": 0.3, "moderate": 0.5, "severe": 0.2}
        ... )
    """
    name: str
    severity_distribution: Dict[str, float] = field(default_factory=lambda: {
        "mild": 0.30,
        "moderate": 0.50,
        "severe": 0.20,
    })
    baseline_biomarker_mean: float = 100.0
    baseline_biomarker_sd: float = 25.0


@dataclass
class VirtualPopulationSpec:
    """
    Specification for virtual population generation.

    Attributes:
        demographics: Demographic distribution
        disease: Disease specification (optional)
        covariate_correlations: Covariate correlations
        seed: Random seed for reproducibility

    Example:
        >>> spec = VirtualPopulationSpec(
        ...     demographics=DemographicSpec(age_mean=50.0),
        ...     seed=42
        ... )
    """
    demographics: DemographicSpec = field(default_factory=DemographicSpec)
    disease: Optional[DiseaseSpec] = None
    covariate_correlations: Dict[Tuple[str, str], float] = field(default_factory=dict)
    seed: int = 12345


@dataclass
class VirtualSubject:
    """
    Virtual subject for trial simulation.

    Attributes:
        id: Subject ID
        age: Age (years)
        weight: Weight (kg)
        sex: Sex ('male' or 'female')
        race: Race/ethnicity
        disease_severity: Disease severity level (optional)
        baseline_biomarker: Baseline biomarker value (optional)
        other: Other covariates

    Example:
        >>> subject = VirtualSubject(
        ...     id=1, age=45.0, weight=75.0, sex='male', race='caucasian'
        ... )
    """
    id: int
    age: float
    weight: float
    sex: str
    race: str
    disease_severity: Optional[str] = None
    baseline_biomarker: Optional[float] = None
    other: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "id": self.id,
            "age": self.age,
            "weight": self.weight,
            "sex": self.sex,
            "race": self.race,
            "disease_severity": self.disease_severity,
            "baseline_biomarker": self.baseline_biomarker,
            "other": self.other,
        }


def _sample_truncated_normal(
    rng: random.Random,
    mean: float,
    sd: float,
    lower: float,
    upper: float,
) -> float:
    """Sample from a truncated normal distribution."""
    for _ in range(1000):
        x = rng.gauss(mean, sd)
        if lower <= x <= upper:
            return x
    return max(lower, min(upper, rng.gauss(mean, sd)))


def _sample_categorical(rng: random.Random, distribution: Dict[str, float]) -> str:
    """Sample from a categorical distribution."""
    u = rng.random()
    cumsum = 0.0
    for category, prob in distribution.items():
        cumsum += prob
        if u <= cumsum:
            return category
    return list(distribution.keys())[-1]


def default_demographic_spec() -> DemographicSpec:
    """
    Create default demographic specification for healthy volunteers.

    Returns:
        DemographicSpec object

    Example:
        >>> spec = default_demographic_spec()
    """
    return DemographicSpec(
        age_mean=35.0,
        age_sd=10.0,
        age_range=(18.0, 55.0),
        weight_mean=75.0,
        weight_sd=12.0,
        weight_range=(50.0, 120.0),
        female_proportion=0.5,
    )


def healthy_volunteer_spec() -> DemographicSpec:
    """
    Create demographic specification for healthy volunteer studies.

    Returns:
        DemographicSpec object

    Example:
        >>> spec = healthy_volunteer_spec()
    """
    return DemographicSpec(
        age_mean=30.0,
        age_sd=8.0,
        age_range=(18.0, 45.0),
        weight_mean=72.0,
        weight_sd=10.0,
        weight_range=(55.0, 100.0),
        female_proportion=0.5,
        race_distribution={
            "caucasian": 0.65,
            "asian": 0.20,
            "black": 0.10,
            "other": 0.05,
        },
    )


def patient_population_spec(
    disease: str,
    severity_weights: Optional[Dict[str, float]] = None,
) -> Tuple[DemographicSpec, DiseaseSpec]:
    """
    Create demographic and disease specification for patient population.

    Args:
        disease: Disease type ('diabetes', 'hepatic', 'renal', 'oncology')
        severity_weights: Custom severity distribution

    Returns:
        Tuple of (DemographicSpec, DiseaseSpec)

    Example:
        >>> demo, dis = patient_population_spec('diabetes')
    """
    if disease == "diabetes":
        demo = DemographicSpec(
            age_mean=55.0, age_sd=12.0, age_range=(30.0, 75.0),
            weight_mean=90.0, weight_sd=18.0, weight_range=(60.0, 150.0),
            female_proportion=0.45,
            race_distribution={"caucasian": 0.60, "asian": 0.15, "black": 0.15, "hispanic": 0.10},
        )
        severity = severity_weights or {"mild": 0.25, "moderate": 0.50, "severe": 0.25}
        dis = DiseaseSpec(
            name="diabetes",
            severity_distribution=severity,
            baseline_biomarker_mean=8.5,  # HbA1c
            baseline_biomarker_sd=1.5,
        )

    elif disease == "hepatic":
        demo = DemographicSpec(
            age_mean=52.0, age_sd=14.0, age_range=(25.0, 75.0),
            weight_mean=80.0, weight_sd=15.0, weight_range=(50.0, 130.0),
            female_proportion=0.40,
        )
        severity = severity_weights or {"mild": 0.40, "moderate": 0.40, "severe": 0.20}
        dis = DiseaseSpec(
            name="hepatic_impairment",
            severity_distribution=severity,
            baseline_biomarker_mean=2.5,
            baseline_biomarker_sd=0.8,
        )

    elif disease == "renal":
        demo = DemographicSpec(
            age_mean=58.0, age_sd=15.0, age_range=(25.0, 80.0),
            weight_mean=78.0, weight_sd=16.0, weight_range=(45.0, 130.0),
            female_proportion=0.45,
            race_distribution={"caucasian": 0.60, "asian": 0.15, "black": 0.20, "other": 0.05},
        )
        severity = severity_weights or {"mild": 0.30, "moderate": 0.40, "severe": 0.20, "esrd": 0.10}
        dis = DiseaseSpec(
            name="renal_impairment",
            severity_distribution=severity,
            baseline_biomarker_mean=45.0,  # eGFR
            baseline_biomarker_sd=20.0,
        )

    elif disease == "oncology":
        demo = DemographicSpec(
            age_mean=62.0, age_sd=12.0, age_range=(25.0, 85.0),
            weight_mean=72.0, weight_sd=15.0, weight_range=(40.0, 120.0),
            female_proportion=0.48,
        )
        severity = severity_weights or {"stage_I": 0.15, "stage_II": 0.25, "stage_III": 0.35, "stage_IV": 0.25}
        dis = DiseaseSpec(
            name="cancer",
            severity_distribution=severity,
            baseline_biomarker_mean=50.0,
            baseline_biomarker_sd=30.0,
        )

    else:
        demo = DemographicSpec(
            age_mean=50.0, age_sd=15.0, age_range=(18.0, 80.0),
            weight_mean=78.0, weight_sd=15.0, weight_range=(45.0, 140.0),
            female_proportion=0.50,
        )
        severity = severity_weights or {"mild": 0.30, "moderate": 0.50, "severe": 0.20}
        dis = DiseaseSpec(
            name=disease,
            severity_distribution=severity,
            baseline_biomarker_mean=100.0,
            baseline_biomarker_sd=25.0,
        )

    return demo, dis


def generate_virtual_population(
    n: int,
    spec: Optional[VirtualPopulationSpec] = None,
    demographics: Optional[DemographicSpec] = None,
    disease: Optional[DiseaseSpec] = None,
    seed: Optional[int] = None,
) -> List[VirtualSubject]:
    """
    Generate a virtual population.

    Args:
        n: Number of subjects
        spec: Complete population specification (overrides other args)
        demographics: Demographic specification
        disease: Disease specification
        seed: Random seed

    Returns:
        List of VirtualSubject objects

    Example:
        >>> pop = generate_virtual_population(100)
        >>> pop = generate_virtual_population(100, seed=42)
        >>> pop = generate_virtual_population(
        ...     100,
        ...     demographics=DemographicSpec(age_mean=50.0)
        ... )
    """
    if spec is not None:
        demographics = spec.demographics
        disease = spec.disease
        seed = spec.seed

    if demographics is None:
        demographics = DemographicSpec()

    if seed is None:
        seed = 12345

    rng = random.Random(seed)
    population = []

    for i in range(n):
        age = _sample_truncated_normal(
            rng, demographics.age_mean, demographics.age_sd,
            demographics.age_range[0], demographics.age_range[1]
        )
        weight = _sample_truncated_normal(
            rng, demographics.weight_mean, demographics.weight_sd,
            demographics.weight_range[0], demographics.weight_range[1]
        )
        sex = "female" if rng.random() < demographics.female_proportion else "male"
        race = _sample_categorical(rng, demographics.race_distribution)

        disease_severity = None
        baseline_biomarker = None

        if disease is not None:
            disease_severity = _sample_categorical(rng, disease.severity_distribution)
            baseline_biomarker = _sample_truncated_normal(
                rng,
                disease.baseline_biomarker_mean,
                disease.baseline_biomarker_sd,
                max(0.0, disease.baseline_biomarker_mean - 3 * disease.baseline_biomarker_sd),
                disease.baseline_biomarker_mean + 3 * disease.baseline_biomarker_sd,
            )

        subject = VirtualSubject(
            id=i + 1,
            age=age,
            weight=weight,
            sex=sex,
            race=race,
            disease_severity=disease_severity,
            baseline_biomarker=baseline_biomarker,
        )
        population.append(subject)

    return population


def summarize_population(population: List[VirtualSubject]) -> Dict[str, Any]:
    """
    Generate summary statistics for a virtual population.

    Args:
        population: List of VirtualSubject objects

    Returns:
        Dictionary of summary statistics

    Example:
        >>> pop = generate_virtual_population(100)
        >>> summary = summarize_population(pop)
        >>> print(f"Mean age: {summary['age_mean']:.1f}")
    """
    n = len(population)
    if n == 0:
        return {"n": 0}

    ages = [s.age for s in population]
    weights = [s.weight for s in population]
    n_female = sum(1 for s in population if s.sex == "female")

    # Race counts
    race_counts: Dict[str, int] = {}
    for s in population:
        race_counts[s.race] = race_counts.get(s.race, 0) + 1

    # Disease severity counts
    severity_counts: Dict[str, int] = {}
    has_disease = any(s.disease_severity is not None for s in population)
    if has_disease:
        for s in population:
            if s.disease_severity is not None:
                severity_counts[s.disease_severity] = severity_counts.get(s.disease_severity, 0) + 1

    age_mean = sum(ages) / n
    weight_mean = sum(weights) / n

    summary = {
        "n": n,
        "age_mean": age_mean,
        "age_sd": math.sqrt(sum((a - age_mean) ** 2 for a in ages) / (n - 1)) if n > 1 else 0.0,
        "age_range": (min(ages), max(ages)),
        "weight_mean": weight_mean,
        "weight_sd": math.sqrt(sum((w - weight_mean) ** 2 for w in weights) / (n - 1)) if n > 1 else 0.0,
        "weight_range": (min(weights), max(weights)),
        "female_proportion": n_female / n,
        "race_distribution": {k: v / n for k, v in race_counts.items()},
    }

    if has_disease:
        biomarkers = [s.baseline_biomarker for s in population if s.baseline_biomarker is not None]
        if biomarkers:
            bm_mean = sum(biomarkers) / len(biomarkers)
            summary["biomarker_mean"] = bm_mean
            summary["biomarker_sd"] = math.sqrt(
                sum((b - bm_mean) ** 2 for b in biomarkers) / (len(biomarkers) - 1)
            ) if len(biomarkers) > 1 else 0.0
        summary["severity_distribution"] = {k: v / n for k, v in severity_counts.items()}

    return summary
