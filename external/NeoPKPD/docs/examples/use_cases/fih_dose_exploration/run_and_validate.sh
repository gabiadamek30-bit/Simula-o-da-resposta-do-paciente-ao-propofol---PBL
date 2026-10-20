#!/usr/bin/env bash
set -euo pipefail

julia docs/examples/use_cases/fih_dose_exploration/run.jl
julia docs/examples/use_cases/fih_dose_exploration/validate.jl
