#!/bin/bash
#
# NeoPKPD Benchmark Suite - Run All Benchmarks
# =============================================
#
# Usage:
#   cd packages/core/benchmarks
#   ./run_all.sh [options]
#
# Options:
#   --julia-only    Only run Julia benchmarks
#   --skip-r        Skip R benchmarks (mrgsolve, nlmixr2)
#   --quick         Reduce n_runs for quick testing
#

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Parse arguments
JULIA_ONLY=false
SKIP_R=false
QUICK=false

for arg in "$@"; do
    case $arg in
        --julia-only)
            JULIA_ONLY=true
            ;;
        --skip-r)
            SKIP_R=true
            ;;
        --quick)
            QUICK=true
            ;;
    esac
done

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORE_DIR="$(dirname "$SCRIPT_DIR")"
REPO_ROOT="$(dirname "$(dirname "$CORE_DIR")")"

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║           NeoPKPD Benchmark Suite - Full Run                    ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Create results directory
mkdir -p "$SCRIPT_DIR/results"
mkdir -p "$SCRIPT_DIR/figures"

# Record system info
echo -e "${YELLOW}System Information:${NC}"
echo "  Date: $(date)"
echo "  Hostname: $(hostname)"
echo "  OS: $(uname -s)"
echo "  Architecture: $(uname -m)"
if command -v julia &> /dev/null; then
    echo "  Julia: $(julia --version)"
fi
if command -v python3 &> /dev/null; then
    echo "  Python: $(python3 --version)"
fi
if command -v R &> /dev/null; then
    echo "  R: $(R --version | head -1)"
fi
echo ""

# ============================================================================
# 1. Julia Benchmarks
# ============================================================================

echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}[1/4] Running NeoPKPD Julia Benchmarks${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if command -v julia &> /dev/null; then
    cd "$CORE_DIR"
    julia --project=. "$SCRIPT_DIR/scripts/benchmark_neopkpd.jl"
    echo -e "${GREEN}✓ Julia benchmarks complete${NC}"
else
    echo -e "${RED}✗ Julia not found, skipping${NC}"
fi

if [ "$JULIA_ONLY" = true ]; then
    echo ""
    echo -e "${YELLOW}--julia-only specified, skipping other benchmarks${NC}"
    exit 0
fi

# ============================================================================
# 2. Python Benchmarks
# ============================================================================

echo ""
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}[2/4] Running NeoPKPD Python Benchmarks${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if command -v python3 &> /dev/null; then
    cd "$REPO_ROOT/packages/python"

    # Check if neopkpd is installed
    if python3 -c "import neopkpd" 2>/dev/null; then
        python3 "$SCRIPT_DIR/scripts/benchmark_neopkpd_python.py"
        echo -e "${GREEN}✓ Python benchmarks complete${NC}"
    else
        echo -e "${YELLOW}! neopkpd not installed, installing...${NC}"
        pip install -e . -q
        python3 "$SCRIPT_DIR/scripts/benchmark_neopkpd_python.py"
        echo -e "${GREEN}✓ Python benchmarks complete${NC}"
    fi
else
    echo -e "${RED}✗ Python not found, skipping${NC}"
fi

# ============================================================================
# 3. R Benchmarks (mrgsolve, nlmixr2)
# ============================================================================

if [ "$SKIP_R" = false ]; then
    echo ""
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}[3/4] Running R Comparison Benchmarks${NC}"
    echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    if command -v Rscript &> /dev/null; then
        cd "$SCRIPT_DIR"

        # mrgsolve
        echo -e "${BLUE}Running mrgsolve benchmarks...${NC}"
        if Rscript -e "library(mrgsolve)" 2>/dev/null; then
            Rscript scripts/benchmark_mrgsolve.R
            echo -e "${GREEN}✓ mrgsolve benchmarks complete${NC}"
        else
            echo -e "${YELLOW}! mrgsolve not installed, skipping${NC}"
        fi

        # nlmixr2
        echo -e "${BLUE}Running nlmixr2 benchmarks...${NC}"
        if Rscript -e "library(nlmixr2)" 2>/dev/null; then
            Rscript scripts/benchmark_nlmixr2.R
            echo -e "${GREEN}✓ nlmixr2 benchmarks complete${NC}"
        else
            echo -e "${YELLOW}! nlmixr2 not installed, skipping${NC}"
        fi
    else
        echo -e "${RED}✗ R not found, skipping${NC}"
    fi
else
    echo ""
    echo -e "${YELLOW}--skip-r specified, skipping R benchmarks${NC}"
fi

# ============================================================================
# 4. Generate Figures
# ============================================================================

echo ""
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}[4/4] Generating Publication Figures${NC}"
echo -e "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

if command -v python3 &> /dev/null; then
    cd "$SCRIPT_DIR"

    # Check for matplotlib
    if python3 -c "import matplotlib" 2>/dev/null; then
        python3 scripts/plot_benchmarks.py
        echo -e "${GREEN}✓ Figures generated${NC}"
    else
        echo -e "${YELLOW}! matplotlib not installed, installing...${NC}"
        pip install matplotlib pandas numpy -q
        python3 scripts/plot_benchmarks.py
        echo -e "${GREEN}✓ Figures generated${NC}"
    fi
else
    echo -e "${RED}✗ Python not found, cannot generate figures${NC}"
fi

# ============================================================================
# Summary
# ============================================================================

echo ""
echo -e "${BLUE}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                    Benchmark Suite Complete                      ║${NC}"
echo -e "${BLUE}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo "Results saved to:"
echo "  $SCRIPT_DIR/results/"
echo ""
echo "Figures saved to:"
echo "  $SCRIPT_DIR/figures/"
echo ""

# List generated files
echo "Generated files:"
ls -la "$SCRIPT_DIR/results/"*.csv 2>/dev/null || echo "  No CSV files found"
echo ""
ls -la "$SCRIPT_DIR/figures/"*.png 2>/dev/null || echo "  No PNG files found"
echo ""
