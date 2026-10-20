# 3+3 Dose Escalation Design

Traditional rule-based dose escalation for Phase I first-in-human studies.

## Design

- **Type**: 3+3 dose escalation
- **Dose levels**: 10, 25, 50, 100, 200 mg
- **Cohort size**: 3 subjects
- **Max DLT rate**: 33%
- **Max subjects**: 30

## Escalation Rules

| DLTs in Cohort | Action |
|----------------|--------|
| 0/3 | Escalate to next dose |
| 1/3 | Expand cohort (+3) |
| 0-1/6 | Escalate to next dose |
| 2+/3 or 2+/6 | Stop, previous dose is MTD |

## Features Demonstrated

- Dose escalation design
- Rule-based decisions
- DLT probability modeling
- MTD determination

## Files

| File | Description |
|------|-------------|
| [python.py](python.py) | Python implementation |
| [julia.jl](julia.jl) | Julia implementation |
| [cli.json](cli.json) | CLI specification |

## Expected Output

```
3+3 Dose Escalation Study
=========================
Dose levels: [10, 25, 50, 100, 200] mg
Starting dose: 10 mg
Max subjects: 30

Cohort 1 (10 mg): 0/3 DLTs -> Escalate
Cohort 2 (25 mg): 0/3 DLTs -> Escalate
Cohort 3 (50 mg): 1/3 DLTs -> Expand
Cohort 4 (50 mg): 0/3 DLTs (1/6 total) -> Escalate
Cohort 5 (100 mg): 2/3 DLTs -> Stop

Recommended MTD: 50 mg
Total subjects: 15
```
