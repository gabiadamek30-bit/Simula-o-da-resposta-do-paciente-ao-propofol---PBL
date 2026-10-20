# Time-Varying Covariates

Population simulation with covariates that change over time (e.g., renal function).

## Model

Creatinine clearance (CrCL) affecting drug clearance:
```
CL(t) = CL_pop × (CrCL(t)/100)^0.5 × exp(η_CL)
```

## Use Cases

- Critically ill patients (changing renal function)
- Chronic therapy (disease progression)
- Drug interactions (enzyme induction over time)
