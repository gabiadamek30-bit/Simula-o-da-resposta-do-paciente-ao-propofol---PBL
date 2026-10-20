#!/usr/bin/env bash
set -euo pipefail

julia docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/check_schema.jl
julia docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/run.jl
julia docs/examples/real_world_validation/datasets/warfarin_nlmixr2data/validate.jl
