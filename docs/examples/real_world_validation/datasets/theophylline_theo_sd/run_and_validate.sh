#!/usr/bin/env bash
set -euo pipefail

julia docs/examples/real_world_validation/datasets/theophylline_theo_sd/run.jl
julia docs/examples/real_world_validation/datasets/theophylline_theo_sd/validate.jl
