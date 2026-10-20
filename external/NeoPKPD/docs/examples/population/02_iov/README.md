# Inter-Occasion Variability (IOV)

Population simulation with both IIV and IOV for multiple dosing occasions.

## Model

```
θᵢⱼ = θ_pop × exp(ηᵢ + κⱼ)

Where:
- ηᵢ ~ N(0, ω²)  (between-subject)
- κⱼ ~ N(0, π²)  (between-occasion)
```

## Parameters

| Parameter | Population | ω (IIV) | π (IOV) |
|-----------|------------|---------|---------|
| CL | 5.0 L/h | 0.3 | 0.15 |
| V | 50.0 L | 0.2 | 0.1 |

## What This Example Shows

1. 50 subjects with 3 dosing occasions each
2. Each subject has consistent IIV across occasions
3. Additional IOV varies by occasion
4. Demonstrates within-subject variability

## Use Cases

- Multiple dose studies
- Crossover bioequivalence
- Chronic therapy monitoring
