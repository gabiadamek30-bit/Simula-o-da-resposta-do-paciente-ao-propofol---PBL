#!/usr/bin/env bash
set -euo pipefail

# Use existing venv if present, otherwise create one
if [ -d "packages/python/.venv" ]; then
    source packages/python/.venv/bin/activate
else
    python3 -m venv packages/python/.venv
    source packages/python/.venv/bin/activate
fi

pip install -e packages/python
pip install pytest

pytest -q packages/python/tests

python packages/python/examples/write_pk_iv_bolus.py
