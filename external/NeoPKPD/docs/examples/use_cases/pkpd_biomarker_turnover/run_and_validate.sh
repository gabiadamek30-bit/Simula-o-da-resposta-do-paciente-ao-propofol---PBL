#!/usr/bin/env bash
set -euo pipefail

julia docs/examples/use_cases/pkpd_biomarker_turnover/run.jl
julia docs/examples/use_cases/pkpd_biomarker_turnover/validate.jl
