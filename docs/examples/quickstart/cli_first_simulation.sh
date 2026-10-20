#!/bin/bash
# NeoPKPD Quickstart - CLI
# Run: ./docs/examples/quickstart/cli_first_simulation.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
OUTPUT_DIR="$SCRIPT_DIR/output"

echo "NeoPKPD Quickstart - CLI"
echo "========================================"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# 1. Run simulation
echo ""
echo "1. Running simulation..."
"$ROOT_DIR/bin/neopkpd" simulate \
    --spec "$SCRIPT_DIR/spec.json" \
    --out "$OUTPUT_DIR/result.json"

echo "   Simulation complete: $OUTPUT_DIR/result.json"

# 2. Compute metrics
echo ""
echo "2. Computing PK metrics..."
"$ROOT_DIR/bin/neopkpd" metrics \
    --artifact "$OUTPUT_DIR/result.json" \
    --metrics cmax,tmax,auc_0_24,auc_0_inf,t_half

# 3. Replay to verify reproducibility
echo ""
echo "3. Verifying reproducibility (replay)..."
"$ROOT_DIR/bin/neopkpd" replay \
    --artifact "$OUTPUT_DIR/result.json" \
    --out "$OUTPUT_DIR/replayed.json"

echo "   Replay successful!"

# 4. Display key results
echo ""
echo "========================================"
echo "RESULTS"
echo "========================================"
echo ""
echo "Output files:"
echo "  - Simulation: $OUTPUT_DIR/result.json"
echo "  - Replayed:   $OUTPUT_DIR/replayed.json"
echo ""
echo "View concentration-time data:"
echo "  jq '.result.observations.conc[:10]' $OUTPUT_DIR/result.json"
echo ""
echo "Quickstart complete!"
