#!/usr/bin/env bash
set -euo pipefail

julia docs/examples/real_world_validation/datasets/theophylline_theo_md/check_schema.jl
julia docs/examples/real_world_validation/datasets/theophylline_theo_md/run.jl
julia docs/examples/real_world_validation/datasets/theophylline_theo_md/validate.jl
