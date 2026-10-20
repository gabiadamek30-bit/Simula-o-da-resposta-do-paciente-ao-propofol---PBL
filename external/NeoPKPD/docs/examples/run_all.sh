#!/usr/bin/env bash
set -euo pipefail

# Ensure output dir exists and is clean
mkdir -p docs/examples/output
rm -f docs/examples/output/*.json

echo "Running Julia docs examples"
julia docs/examples/julia/01_pk_iv_bolus.jl

echo "Running Python docs examples"
python3 -m pip install -e packages/python
python3 docs/examples/python/01_replay_golden.py

echo "Replaying generated doc artifacts"
./packages/cli/bin/neopkpd replay --artifact docs/examples/output/01_pk_iv_bolus.json

echo "Validating doc outputs metadata"
julia docs/examples/validate_outputs.jl
