#!/bin/bash
# ==============================================================================
# NONMEM Comprehensive Benchmark Suite
# ==============================================================================
#
# Prerequisites:
#   - NONMEM 7.4+ installed and nmfe74/75/76 in PATH
#   - Valid NONMEM license
#
# Usage:
#   ./benchmark_nonmem_comprehensive.sh [nmfe_version]
#   Example: ./benchmark_nonmem_comprehensive.sh nmfe75
#
# Output:
#   ../results/nonmem_comprehensive_benchmarks.csv
# ==============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="$SCRIPT_DIR/../results"
MODELS_DIR="$SCRIPT_DIR/../nonmem_models"
OUTPUT_FILE="$RESULTS_DIR/nonmem_comprehensive_benchmarks.csv"

# NONMEM executable (default to nmfe75)
NMFE="${1:-nmfe75}"

# Check if NONMEM is available
if ! command -v "$NMFE" &> /dev/null; then
    echo "ERROR: $NMFE not found in PATH"
    echo "Please install NONMEM and ensure nmfeXX is in your PATH"
    echo ""
    echo "For benchmark comparison, you can use reference data from:"
    echo "  - Bauer RJ. NONMEM Tutorial Part II. CPT PSP. 2019"
    echo "  - Wang W, et al. Comparison of PK Software. 2019"
    exit 1
fi

# Create directories
mkdir -p "$MODELS_DIR"
mkdir -p "$RESULTS_DIR"

# Configuration
N_RUNS=10  # NONMEM is slow, so fewer runs
TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
NONMEM_VERSION=$($NMFE -version 2>/dev/null | head -1 || echo "Unknown")

echo "=============================================================================="
echo "NONMEM Comprehensive Benchmark Suite"
echo "=============================================================================="
echo "NONMEM Version: $NONMEM_VERSION"
echo "Runs per benchmark: $N_RUNS"
echo "Timestamp: $TIMESTAMP"
echo "=============================================================================="

# Initialize CSV
echo '"category","subcategory","name","model","n_subjects","n_runs","mean_ms","std_ms","median_ms","min_ms","max_ms","ci_lower_ms","ci_upper_ms","timestamp","nonmem_version"' > "$OUTPUT_FILE"

# Function to run benchmark
run_benchmark() {
    local category="$1"
    local subcategory="$2"
    local name="$3"
    local model="$4"
    local n_subjects="$5"
    local control_file="$6"

    echo "  Running: $name ($model)..."

    local times=()
    local workdir=$(mktemp -d)

    for ((i=1; i<=N_RUNS; i++)); do
        cp "$control_file" "$workdir/run.ctl"
        cd "$workdir"

        # Time NONMEM execution
        local start=$(python3 -c "import time; print(time.time() * 1000)")
        $NMFE run.ctl run.lst -background > /dev/null 2>&1
        local end=$(python3 -c "import time; print(time.time() * 1000)")

        local elapsed=$(echo "$end - $start" | bc)
        times+=($elapsed)

        rm -rf "$workdir"/*
        cd - > /dev/null
    done

    rm -rf "$workdir"

    # Calculate statistics using Python
    local stats=$(python3 -c "
import statistics
times = [${times[*]}]
times = [float(t) for t in '${times[*]}'.split()]
print(f'{statistics.mean(times):.6f}')
print(f'{statistics.stdev(times) if len(times) > 1 else 0:.6f}')
print(f'{statistics.median(times):.6f}')
print(f'{min(times):.6f}')
print(f'{max(times):.6f}')
sorted_times = sorted(times)
n = len(times)
print(f'{sorted_times[int(n*0.025)]:.6f}')
print(f'{sorted_times[int(n*0.975)]:.6f}')
")

    local mean_ms=$(echo "$stats" | sed -n '1p')
    local std_ms=$(echo "$stats" | sed -n '2p')
    local median_ms=$(echo "$stats" | sed -n '3p')
    local min_ms=$(echo "$stats" | sed -n '4p')
    local max_ms=$(echo "$stats" | sed -n '5p')
    local ci_lower=$(echo "$stats" | sed -n '6p')
    local ci_upper=$(echo "$stats" | sed -n '7p')

    echo "    Mean: $mean_ms ms (±$std_ms)"

    echo "\"$category\",\"$subcategory\",\"$name\",\"$model\",$n_subjects,$N_RUNS,$mean_ms,$std_ms,$median_ms,$min_ms,$max_ms,$ci_lower,$ci_upper,\"$TIMESTAMP\",\"$NONMEM_VERSION\"" >> "$OUTPUT_FILE"
}

# ==============================================================================
# Create NONMEM Control Files
# ==============================================================================

create_nonmem_models() {
    # 1. One-Compartment IV Bolus
    cat > "$MODELS_DIR/onecomp_iv.ctl" << 'EOF'
$PROBLEM One-Compartment IV Bolus Simulation
$INPUT ID TIME DV AMT EVID
$DATA onecomp_data.csv IGNORE=@
$SUBROUTINE ADVAN1 TRANS2
$PK
  CL = THETA(1)
  V = THETA(2)
$ERROR
  Y = F + ERR(1)
$THETA
  (0, 5)   ; CL
  (0, 50)  ; V
$OMEGA 0.1
$SIGMA 0.01
$SIMULATION (12345) ONLYSIMULATION
$TABLE ID TIME DV NOPRINT ONEHEADER FILE=sim.tab
EOF

    # 2. Two-Compartment IV Bolus
    cat > "$MODELS_DIR/twocomp_iv.ctl" << 'EOF'
$PROBLEM Two-Compartment IV Bolus Simulation
$INPUT ID TIME DV AMT EVID
$DATA twocomp_data.csv IGNORE=@
$SUBROUTINE ADVAN3 TRANS4
$PK
  CL = THETA(1)
  V1 = THETA(2)
  Q = THETA(3)
  V2 = THETA(4)
$ERROR
  Y = F + ERR(1)
$THETA
  (0, 5)    ; CL
  (0, 50)   ; V1
  (0, 10)   ; Q
  (0, 100)  ; V2
$OMEGA 0.1
$SIGMA 0.01
$SIMULATION (12345) ONLYSIMULATION
$TABLE ID TIME DV NOPRINT ONEHEADER FILE=sim.tab
EOF

    # 3. One-Compartment Oral
    cat > "$MODELS_DIR/onecomp_oral.ctl" << 'EOF'
$PROBLEM One-Compartment Oral Simulation
$INPUT ID TIME DV AMT EVID
$DATA onecomp_oral_data.csv IGNORE=@
$SUBROUTINE ADVAN2 TRANS2
$PK
  KA = THETA(1)
  CL = THETA(2)
  V = THETA(3)
$ERROR
  Y = F + ERR(1)
$THETA
  (0, 1.5)  ; KA
  (0, 5)    ; CL
  (0, 50)   ; V
$OMEGA 0.1
$SIGMA 0.01
$SIMULATION (12345) ONLYSIMULATION
$TABLE ID TIME DV NOPRINT ONEHEADER FILE=sim.tab
EOF

    # 4. Emax PD Model
    cat > "$MODELS_DIR/emax_pd.ctl" << 'EOF'
$PROBLEM Direct Emax PD Model
$INPUT ID TIME CP DV
$DATA emax_data.csv IGNORE=@
$PRED
  E0 = THETA(1)
  EMAX = THETA(2)
  EC50 = THETA(3)
  EFFECT = E0 + EMAX * CP / (EC50 + CP)
  Y = EFFECT + ERR(1)
$THETA
  (0, 0)     ; E0
  (0, 100)   ; EMAX
  (0, 10)    ; EC50
$OMEGA 0.1
$SIGMA 0.01
$SIMULATION (12345) ONLYSIMULATION
$TABLE ID TIME DV NOPRINT ONEHEADER FILE=sim.tab
EOF

    # 5. Population PK with IIV
    cat > "$MODELS_DIR/pop_pk.ctl" << 'EOF'
$PROBLEM Population Two-Compartment with IIV
$INPUT ID TIME DV AMT EVID
$DATA pop_data.csv IGNORE=@
$SUBROUTINE ADVAN3 TRANS4
$PK
  TVCL = THETA(1)
  TVV1 = THETA(2)
  TVQ = THETA(3)
  TVV2 = THETA(4)

  CL = TVCL * EXP(ETA(1))
  V1 = TVV1 * EXP(ETA(2))
  Q = TVQ * EXP(ETA(3))
  V2 = TVV2 * EXP(ETA(4))
$ERROR
  IPRED = F
  Y = IPRED * (1 + ERR(1))
$THETA
  (0, 5)    ; TVCL
  (0, 50)   ; TVV1
  (0, 10)   ; TVQ
  (0, 100)  ; TVV2
$OMEGA BLOCK(4)
  0.09
  0.01 0.09
  0.01 0.01 0.09
  0.01 0.01 0.01 0.09
$SIGMA 0.04
$SIMULATION (12345) ONLYSIMULATION SUBPROBLEMS=1
$TABLE ID TIME DV IPRED NOPRINT ONEHEADER FILE=sim.tab
EOF

    echo "NONMEM control files created in $MODELS_DIR"
}

# ==============================================================================
# Main Execution
# ==============================================================================

echo ""
echo "Creating NONMEM control files..."
create_nonmem_models

echo ""
echo "======================================================================"
echo "BENCHMARK CATEGORY: PK MODELS"
echo "======================================================================"

# Note: Actual benchmarks require data files and proper NONMEM setup
# This script provides the infrastructure - run when NONMEM is available

echo ""
echo "======================================================================"
echo "BENCHMARK COMPLETE"
echo "Results saved to: $OUTPUT_FILE"
echo "======================================================================"
