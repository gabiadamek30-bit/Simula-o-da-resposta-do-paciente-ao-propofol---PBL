# Virtual Population

Comprehensive guide for generating realistic virtual patient populations in NeoPKPD Python.

---

## Overview

Virtual population generation creates demographically realistic subjects for trial simulation.

```python
from neopkpd import trial

# Generate healthy volunteers
pop = trial.generate_virtual_population(
    n=100,
    spec=trial.healthy_volunteer_spec(),
    seed=42
)

# Access demographics
print(f"Mean age: {pop.mean_age:.1f} years")
print(f"Mean weight: {pop.mean_weight:.1f} kg")
```

---

## DemographicSpec Class

### Class Definition

```python
@dataclass
class DemographicSpec:
    """Specification for virtual population demographics."""

    # Age distribution
    age_mean: float = 35.0
    age_sd: float = 10.0
    age_min: float = 18.0
    age_max: float = 65.0

    # Weight distribution
    weight_mean: float = 75.0
    weight_sd: float = 15.0
    weight_min: float = 45.0
    weight_max: float = 150.0

    # Height distribution
    height_mean: float = 170.0
    height_sd: float = 10.0
    height_min: float = 145.0
    height_max: float = 200.0

    # Sex distribution
    female_fraction: float = 0.50

    # Race distribution
    race_distribution: dict[str, float] | None = None

    # Renal function
    egfr_distribution: dict[str, float] | None = None

    # Hepatic function
    child_pugh_distribution: dict[str, float] | None = None

    # Custom covariates
    covariates: dict[str, CovariateSpec] | None = None
```

### Creating Custom Specs

```python
from neopkpd.trial import DemographicSpec

# Custom population
spec = DemographicSpec(
    age_mean=55.0,
    age_sd=12.0,
    age_min=18.0,
    age_max=85.0,
    weight_mean=85.0,
    weight_sd=18.0,
    weight_min=50.0,
    weight_max=150.0,
    female_fraction=0.45,
    race_distribution={
        "white": 0.60,
        "black": 0.20,
        "asian": 0.15,
        "other": 0.05
    }
)

pop = trial.generate_virtual_population(n=200, spec=spec, seed=42)
```

---

## Pre-built Population Specs

### Healthy Volunteers

```python
# Phase I typical population
spec = trial.healthy_volunteer_spec()

# Parameters:
# Age: 30 ± 8 years (18-45)
# Weight: 72 ± 10 kg (55-90)
# Height: 170 ± 9 cm
# 50% female
# Normal renal/hepatic function
```

### Patient Populations

```python
# General patient population
spec = trial.patient_population_spec()

# Diabetic patients
spec = trial.patient_population_spec(disease="diabetes")
# Includes: HbA1c distribution, older age

# Cardiac patients
spec = trial.patient_population_spec(disease="cardiac")
# Includes: Ejection fraction, cardiac medications

# Oncology patients
spec = trial.patient_population_spec(disease="oncology")
# Includes: Performance status, prior treatments
```

### Special Populations

```python
# Elderly patients
spec = trial.elderly_patient_spec()
# Age: 72 ± 8 years (65-90)
# Reduced renal function

# Pediatric population
spec = trial.pediatric_spec(age_group="child")  # 6-12 years
spec = trial.pediatric_spec(age_group="adolescent")  # 12-18 years
spec = trial.pediatric_spec(age_group="infant")  # 1-6 years

# Renal impairment
spec = trial.renal_impairment_spec(severity="moderate")
# Includes: eGFR 30-59 mL/min/1.73m²

# Hepatic impairment
spec = trial.hepatic_impairment_spec(child_pugh="B")
# Includes: Child-Pugh score B
```

---

## Generating Populations

### Basic Generation

```python
# Generate population
pop = trial.generate_virtual_population(
    n=100,
    spec=spec,
    seed=42
)

# Access as DataFrame
df = pop.to_dataframe()
print(df.head())
```

### VirtualPopulation Class

```python
@dataclass
class VirtualPopulation:
    """Container for virtual population data."""

    n_subjects: int
    subjects: list[VirtualSubject]

    # Summary statistics
    @property
    def mean_age(self) -> float: ...
    @property
    def mean_weight(self) -> float: ...
    @property
    def female_fraction(self) -> float: ...

    def to_dataframe(self) -> pd.DataFrame: ...
    def subset(self, indices: list[int]) -> VirtualPopulation: ...
    def filter(self, criterion: Callable) -> VirtualPopulation: ...
```

### VirtualSubject Class

```python
@dataclass
class VirtualSubject:
    """Individual virtual subject."""

    id: int
    age: float
    weight: float
    height: float
    bmi: float
    sex: str
    race: str
    egfr: float | None = None
    child_pugh: str | None = None
    covariates: dict[str, Any] | None = None
```

---

## Covariate Specification

### Standard Covariates

```python
from neopkpd.trial import CovariateSpec

# Continuous covariate
creatinine = CovariateSpec(
    name="creatinine",
    type="continuous",
    distribution="normal",
    mean=1.0,
    sd=0.3,
    min=0.5,
    max=2.5,
    unit="mg/dL"
)

# Categorical covariate
smoking = CovariateSpec(
    name="smoking_status",
    type="categorical",
    categories=["never", "former", "current"],
    probabilities=[0.50, 0.30, 0.20]
)
```

### Adding Covariates to Spec

```python
# Add custom covariates
spec = DemographicSpec(
    age_mean=55.0,
    age_sd=12.0,
    covariates={
        "creatinine": CovariateSpec(
            type="continuous",
            distribution="lognormal",
            mean=1.0,
            sd=0.3
        ),
        "albumin": CovariateSpec(
            type="continuous",
            distribution="normal",
            mean=4.0,
            sd=0.5,
            min=2.5,
            max=5.5
        ),
        "genotype": CovariateSpec(
            type="categorical",
            categories=["PM", "IM", "EM", "UM"],
            probabilities=[0.05, 0.15, 0.70, 0.10]
        )
    }
)
```

### Correlated Covariates

```python
# Define correlations
spec = DemographicSpec(
    age_mean=55.0,
    age_sd=12.0,
    weight_mean=80.0,
    weight_sd=15.0,
    covariate_correlations={
        ("age", "egfr"): -0.4,      # Age negatively correlated with eGFR
        ("weight", "bmi"): 0.8,      # Weight positively correlated with BMI
        ("age", "creatinine"): 0.3   # Age positively correlated with creatinine
    }
)
```

---

## Renal Function

### eGFR Distribution

```python
# Normal renal function
spec = DemographicSpec(
    egfr_distribution={
        "normal": 1.0,     # eGFR >= 90
        "mild": 0.0,       # eGFR 60-89
        "moderate": 0.0,   # eGFR 30-59
        "severe": 0.0      # eGFR 15-29
    }
)

# Mixed population
spec = DemographicSpec(
    egfr_distribution={
        "normal": 0.40,
        "mild": 0.30,
        "moderate": 0.20,
        "severe": 0.10
    }
)
```

### eGFR Calculation

```python
# eGFR calculated using CKD-EPI equation
for subject in pop.subjects:
    print(f"Subject {subject.id}: eGFR = {subject.egfr:.1f} mL/min/1.73m²")
```

### Renal Impairment Study

```python
# Generate matched groups
def generate_renal_study_population(n_per_group: int, seed: int = 42):
    populations = {}

    for severity in ["normal", "mild", "moderate", "severe"]:
        spec = trial.renal_impairment_spec(severity=severity)
        populations[severity] = trial.generate_virtual_population(
            n=n_per_group,
            spec=spec,
            seed=seed
        )

    return populations

pops = generate_renal_study_population(8)
```

---

## Hepatic Function

### Child-Pugh Distribution

```python
# Hepatic impairment study
spec = DemographicSpec(
    child_pugh_distribution={
        "normal": 0.25,    # No impairment
        "A": 0.25,         # Mild (5-6 points)
        "B": 0.25,         # Moderate (7-9 points)
        "C": 0.25          # Severe (10-15 points)
    }
)
```

### Child-Pugh Calculation

```python
# Child-Pugh score components
@dataclass
class ChildPughComponents:
    bilirubin: float      # mg/dL
    albumin: float        # g/dL
    inr: float            # INR
    ascites: str          # "none", "mild", "moderate"
    encephalopathy: str   # "none", "grade1-2", "grade3-4"

    def calculate_score(self) -> int:
        """Calculate total Child-Pugh score."""
        ...

    def classify(self) -> str:
        """Return A, B, or C classification."""
        ...
```

---

## Population Summary

### Summary Statistics

```python
# Generate summary
summary = trial.summarize_population(pop)

print("=== Population Summary ===")
print(f"N: {summary['n']}")
print(f"Age: {summary['age']['mean']:.1f} ± {summary['age']['sd']:.1f} years")
print(f"Weight: {summary['weight']['mean']:.1f} ± {summary['weight']['sd']:.1f} kg")
print(f"BMI: {summary['bmi']['mean']:.1f} ± {summary['bmi']['sd']:.1f} kg/m²")
print(f"Female: {summary['female_fraction']*100:.0f}%")

# Race breakdown
print("\nRace Distribution:")
for race, frac in summary['race_distribution'].items():
    print(f"  {race}: {frac*100:.0f}%")
```

### Detailed Statistics

```python
# Get detailed stats
stats = trial.population_statistics(pop)

# Continuous variables
for var in ["age", "weight", "height", "bmi"]:
    s = stats[var]
    print(f"{var}: {s['mean']:.1f} [{s['min']:.1f}-{s['max']:.1f}]")

# Categorical variables
for var in ["sex", "race"]:
    print(f"\n{var}:")
    for cat, count in stats[var].items():
        print(f"  {cat}: {count}")
```

---

## Population Filtering

### Filter by Criteria

```python
# Filter by age
elderly = pop.filter(lambda s: s.age >= 65)
print(f"Elderly subjects: {elderly.n_subjects}")

# Filter by multiple criteria
subpop = pop.filter(lambda s: s.age >= 50 and s.sex == "female")

# Filter by renal function
moderate_ri = pop.filter(lambda s: 30 <= s.egfr < 60)
```

### Stratified Sampling

```python
# Sample stratified by sex
stratified = trial.stratified_sample(
    population=pop,
    strata_var="sex",
    n_per_stratum={"male": 25, "female": 25},
    seed=42
)

# Sample stratified by age group
stratified = trial.stratified_sample(
    population=pop,
    strata_var="age_group",
    strata_definition={
        "young": lambda s: s.age < 40,
        "middle": lambda s: 40 <= s.age < 60,
        "elderly": lambda s: s.age >= 60
    },
    n_per_stratum={"young": 20, "middle": 30, "elderly": 20}
)
```

---

## Randomization

### Assign to Arms

```python
# Simple randomization
assignments = trial.randomize(
    population=pop,
    arms=["Placebo", "Treatment"],
    ratio=[1, 1],
    seed=42
)

# Stratified randomization
assignments = trial.stratified_randomize(
    population=pop,
    arms=["Placebo", "Low", "High"],
    ratio=[1, 1, 1],
    strata=["sex", "age_group"],
    block_size=6,
    seed=42
)

# Check balance
balance = trial.check_randomization_balance(assignments)
print(f"Balance metrics: {balance}")
```

---

## Complete Example

```python
from neopkpd import trial

# =========================================
# Phase III Population Generation
# =========================================

# 1. Define target population characteristics
spec = trial.DemographicSpec(
    # Age: typical patient population
    age_mean=58.0,
    age_sd=14.0,
    age_min=18.0,
    age_max=85.0,

    # Weight: slightly overweight
    weight_mean=82.0,
    weight_sd=18.0,
    weight_min=45.0,
    weight_max=160.0,

    # Sex distribution
    female_fraction=0.48,

    # Race distribution (US-based trial)
    race_distribution={
        "white": 0.65,
        "black": 0.18,
        "asian": 0.08,
        "hispanic": 0.06,
        "other": 0.03
    },

    # Include some renal impairment
    egfr_distribution={
        "normal": 0.50,
        "mild": 0.35,
        "moderate": 0.15,
        "severe": 0.0  # Excluded
    },

    # Custom covariates
    covariates={
        "baseline_hba1c": trial.CovariateSpec(
            type="continuous",
            distribution="normal",
            mean=8.2,
            sd=1.2,
            min=7.0,
            max=12.0
        ),
        "diabetes_duration": trial.CovariateSpec(
            type="continuous",
            distribution="lognormal",
            mean=8.0,
            sd=5.0,
            min=1.0,
            max=30.0
        ),
        "prior_therapy": trial.CovariateSpec(
            type="categorical",
            categories=["metformin_only", "dual_therapy", "insulin"],
            probabilities=[0.40, 0.45, 0.15]
        )
    }
)

# 2. Generate population
print("Generating virtual population...")
pop = trial.generate_virtual_population(
    n=300,
    spec=spec,
    seed=42
)

# 3. Display summary
summary = trial.summarize_population(pop)

print("\n=== Population Summary ===")
print(f"Total subjects: {summary['n']}")
print(f"\nDemographics:")
print(f"  Age: {summary['age']['mean']:.1f} ± {summary['age']['sd']:.1f} years")
print(f"  Weight: {summary['weight']['mean']:.1f} ± {summary['weight']['sd']:.1f} kg")
print(f"  BMI: {summary['bmi']['mean']:.1f} ± {summary['bmi']['sd']:.1f} kg/m²")
print(f"  Female: {summary['female_fraction']*100:.0f}%")

print(f"\nRace Distribution:")
for race, frac in summary['race_distribution'].items():
    print(f"  {race}: {frac*100:.0f}%")

print(f"\nRenal Function:")
for category, frac in summary['egfr_distribution'].items():
    print(f"  {category}: {frac*100:.0f}%")

print(f"\nBaseline Characteristics:")
print(f"  HbA1c: {summary['covariates']['baseline_hba1c']['mean']:.1f}%")
print(f"  Diabetes duration: {summary['covariates']['diabetes_duration']['mean']:.1f} years")

# 4. Randomize to treatment arms
assignments = trial.stratified_randomize(
    population=pop,
    arms=["Placebo", "Low Dose", "High Dose"],
    ratio=[1, 1, 1],
    strata=["sex"],
    block_size=6,
    seed=42
)

print(f"\nRandomization:")
for arm, subjects in assignments.items():
    n = len(subjects)
    print(f"  {arm}: {n} subjects")

# 5. Check covariate balance
balance = trial.check_covariate_balance(
    population=pop,
    assignments=assignments,
    covariates=["age", "weight", "baseline_hba1c"]
)

print(f"\nCovariate Balance (p-values):")
for cov, p in balance.items():
    status = "✓" if p > 0.05 else "⚠"
    print(f"  {cov}: p={p:.3f} {status}")
```

---

## See Also

- [Study Designs](designs.md) - Trial design types
- [Dosing Regimens](dosing.md) - Dosing schedules
- [Power Analysis](power.md) - Sample size calculation
- [Julia Population](../../julia/population/index.md) - Julia interface
