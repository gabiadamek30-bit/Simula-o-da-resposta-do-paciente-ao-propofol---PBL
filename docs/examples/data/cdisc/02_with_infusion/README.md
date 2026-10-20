# CDISC Infusion Data Import

Import pharmacokinetic data with IV infusion dosing using EXDUR field.

## Key Fields

### EX Domain for Infusion
- EXDOSE: Dose amount
- EXDOSU: Dose units
- EXDUR: Infusion duration (e.g., "PT1H" for 1 hour)
- EXROUTE: INTRAVENOUS

## Files

| File | Description |
|------|-------------|
| ex.csv | EX domain with EXDUR |
| load.py | Python loader |

Note: This example reuses PC data from the basic example.

## Infusion Handling

```python
# EXDUR format: ISO 8601 duration
# PT1H = 1 hour infusion
# PT30M = 30 minute infusion
# PT2H30M = 2.5 hour infusion

data = load_cdisc(pc="pc.csv", ex="ex.csv")

# Infusion detected from EXDUR
for dose in data.doses:
    print(f"Dose: {dose.amount} mg, Duration: {dose.duration} h")
```
