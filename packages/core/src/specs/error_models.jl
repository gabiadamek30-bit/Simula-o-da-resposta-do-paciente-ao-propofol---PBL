# Residual Error Model Specifications
# Provides additive, proportional, combined, and exponential error models

export ResidualErrorKind, AdditiveError, ProportionalError, CombinedError, ExponentialError
export ResidualErrorParams, AdditiveErrorParams, ProportionalErrorParams
export CombinedErrorParams, ExponentialErrorParams
export ResidualErrorSpec

# -------------------------
# Error Model Kinds
# -------------------------

abstract type ResidualErrorKind end

"""
Additive (constant) error model.

Y = F + ε
ε ~ Normal(0, σ²)

Where:
- Y is the observed value
- F is the predicted value
- σ is the additive standard deviation

The error variance is constant regardless of the predicted value.
Common for assay errors near detection limits.
"""
struct AdditiveError <: ResidualErrorKind end

"""
Proportional (constant CV) error model.

Y = F * (1 + ε)
ε ~ Normal(0, σ²)

Equivalently: Y = F + F*ε, so SD(Y|F) = σ*F

Where:
- Y is the observed value
- F is the predicted value
- σ is the proportional coefficient of variation (CV)

The error standard deviation is proportional to the prediction.
Common for bioanalytical assays.
"""
struct ProportionalError <: ResidualErrorKind end

"""
Combined additive and proportional error model.

Y = F + √(σ_add² + (σ_prop * F)²) * ε
ε ~ Normal(0, 1)

Where:
- σ_add is the additive standard deviation component
- σ_prop is the proportional coefficient of variation

This is the most flexible standard error model, recommended for
most PK/PD applications. It handles both low concentration errors
(dominated by additive) and high concentration errors (dominated
by proportional).
"""
struct CombinedError <: ResidualErrorKind end

"""
Exponential (log-additive) error model.

log(Y) = log(F) + ε
ε ~ Normal(0, σ²)

Equivalently: Y = F * exp(ε)

Where:
- σ is the standard deviation on the log scale

This model ensures predictions are always positive and is
equivalent to assuming Y has a log-normal distribution.
Commonly used for concentration data.
"""
struct ExponentialError <: ResidualErrorKind end

# -------------------------
# Error Model Parameters
# -------------------------

struct AdditiveErrorParams
    sigma::Float64  # Additive SD

    function AdditiveErrorParams(sigma::Float64)
        sigma > 0 || error("sigma must be positive, got $sigma")
        new(sigma)
    end
end

struct ProportionalErrorParams
    sigma::Float64  # Proportional CV

    function ProportionalErrorParams(sigma::Float64)
        sigma > 0 || error("sigma must be positive, got $sigma")
        new(sigma)
    end
end

struct CombinedErrorParams
    sigma_add::Float64   # Additive SD component
    sigma_prop::Float64  # Proportional CV component

    function CombinedErrorParams(sigma_add::Float64, sigma_prop::Float64)
        sigma_add >= 0 || error("sigma_add must be non-negative, got $sigma_add")
        sigma_prop >= 0 || error("sigma_prop must be non-negative, got $sigma_prop")
        (sigma_add > 0 || sigma_prop > 0) || error("At least one sigma must be positive")
        new(sigma_add, sigma_prop)
    end
end

struct ExponentialErrorParams
    sigma::Float64  # Log-scale SD

    function ExponentialErrorParams(sigma::Float64)
        sigma > 0 || error("sigma must be positive, got $sigma")
        new(sigma)
    end
end

# -------------------------
# Residual Error Specification
# -------------------------

"""
Complete specification for residual error model.

Fields:
- kind: Type of error model (AdditiveError, ProportionalError, etc.)
- params: Parameters for the error model
- observation: Which observation to apply error to (e.g., :conc, :effect)
- seed: RNG seed for reproducibility

Examples:
```julia
# Additive error with SD = 0.1
error_spec = ResidualErrorSpec(AdditiveError(), AdditiveErrorParams(0.1), :conc, 12345)

# Combined error (typical for PK)
error_spec = ResidualErrorSpec(
    CombinedError(),
    CombinedErrorParams(0.1, 0.1),  # add=0.1, prop=10% CV
    :conc,
    UInt64(42)
)
```
"""
struct ResidualErrorSpec{K<:ResidualErrorKind,P}
    kind::K
    params::P
    observation::Symbol
    seed::UInt64

    function ResidualErrorSpec(kind::K, params::P, observation::Symbol, seed::Integer) where {K<:ResidualErrorKind,P}
        new{K,P}(kind, params, observation, UInt64(seed))
    end
end
