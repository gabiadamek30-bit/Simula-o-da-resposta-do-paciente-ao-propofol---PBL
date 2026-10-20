#!/usr/bin/env bash
set -euo pipefail

# First instantiate the core package
julia -e 'using Pkg; Pkg.activate("packages/core"); Pkg.instantiate()'

# Develop the local NeoPKPD package and instantiate CLI dependencies
julia -e '
using Pkg
Pkg.activate("packages/cli")
Pkg.develop(path="packages/core")
Pkg.instantiate()
'

./packages/cli/bin/neopkpd version
./packages/cli/bin/neopkpd replay --artifact validation/golden/pk_iv_bolus.json
./packages/cli/bin/neopkpd validate-golden
