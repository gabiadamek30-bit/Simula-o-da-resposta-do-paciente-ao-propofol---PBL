#!/usr/bin/env Rscript
# ==============================================================================
# Monolix Comprehensive Benchmark Suite
# ==============================================================================
#
# Prerequisites:
#   - Monolix installed (2021R1 or later)
#   - lixoftConnectors R package installed
#   - Valid Monolix license
#
# Installation:
#   install.packages("lixoftConnectors", repos = "https://monolix.lixoft.com/R-packages/")
#   lixoftConnectors::initializeLixoftConnectors(software = "monolix")
#
# Usage:
#   Rscript benchmark_monolix_comprehensive.R
#
# Output:
#   ../results/monolix_comprehensive_benchmarks.csv
# ==============================================================================

library(microbenchmark)

# Check for lixoftConnectors
if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
    cat("ERROR: lixoftConnectors package not found\n")
    cat("\nTo install:\n")
    cat("  install.packages('lixoftConnectors', repos = 'https://monolix.lixoft.com/R-packages/')\n")
    cat("\nThen initialize:\n")
    cat("  lixoftConnectors::initializeLixoftConnectors(software = 'monolix')\n")
    cat("\nFor benchmark comparison without Monolix, use reference data from:\n")
    cat("  - Monolix documentation: https://monolix.lixoft.com/\n")
    cat("  - Lavielle M, et al. Monolix Performance. PAGE 2020\n")
    quit(status = 1)
}

library(lixoftConnectors)

# Try to initialize Monolix
tryCatch({
    initializeLixoftConnectors(software = "monolix")
}, error = function(e) {
    cat("ERROR: Could not initialize Monolix\n")
    cat("Please ensure Monolix is installed and licensed.\n")
    quit(status = 1)
})

# ==============================================================================
# Configuration
# ==============================================================================

SCRIPT_DIR <- dirname(sys.frame(1)$ofile)
if (is.null(SCRIPT_DIR)) SCRIPT_DIR <- "."
RESULTS_DIR <- file.path(SCRIPT_DIR, "..", "results")
MODELS_DIR <- file.path(SCRIPT_DIR, "..", "monolix_models")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(MODELS_DIR, showWarnings = FALSE, recursive = TRUE)

CONFIG <- list(
    n_runs = 20,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    monolix_version = getMonolixVersion()
)

cat("==============================================================================\n")
cat("Monolix Comprehensive Benchmark Suite\n")
cat("==============================================================================\n")
cat("Monolix Version:", CONFIG$monolix_version, "\n")
cat("Runs per benchmark:", CONFIG$n_runs, "\n")
cat("Timestamp:", CONFIG$timestamp, "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# Benchmark Infrastructure
# ==============================================================================

results <- data.frame()

run_benchmark <- function(f, name, category, subcategory = "", model, n_subjects = 1) {
    cat(sprintf("  Running: %s (%s)...\n", name, model))

    times <- microbenchmark(f(), times = CONFIG$n_runs, unit = "ms")$time / 1e6  # Convert to ms

    result <- data.frame(
        category = category,
        subcategory = subcategory,
        name = name,
        model = model,
        n_subjects = n_subjects,
        n_runs = CONFIG$n_runs,
        mean_ms = mean(times),
        std_ms = sd(times),
        median_ms = median(times),
        min_ms = min(times),
        max_ms = max(times),
        ci_lower_ms = quantile(times, 0.025),
        ci_upper_ms = quantile(times, 0.975),
        timestamp = CONFIG$timestamp,
        monolix_version = CONFIG$monolix_version
    )

    cat(sprintf("    Mean: %.3f ms (±%.3f)\n", result$mean_ms, result$std_ms))
    return(result)
}

# ==============================================================================
# Create Monolix Model Files
# ==============================================================================

create_monolix_models <- function() {
    # One-Compartment IV Bolus (Mlxtran)
    writeLines('
[LONGITUDINAL]
input = {V, Cl}

EQUATION:
; PK model
Cc = pkmodel(V, Cl)

OUTPUT:
output = Cc
', file.path(MODELS_DIR, "onecomp_iv.txt"))

    # Two-Compartment IV Bolus
    writeLines('
[LONGITUDINAL]
input = {V1, V2, Q, Cl}

EQUATION:
; Two-compartment PK model
Cc = pkmodel(V1, V2, Q, Cl)

OUTPUT:
output = Cc
', file.path(MODELS_DIR, "twocomp_iv.txt"))

    # Emax PD Model
    writeLines('
[LONGITUDINAL]
input = {E0, Emax, EC50}

EQUATION:
; Direct Emax model
E = E0 + Emax * Cc / (EC50 + Cc)

OUTPUT:
output = E
', file.path(MODELS_DIR, "emax_pd.txt"))

    # Indirect Response Model (IRM1)
    writeLines('
[LONGITUDINAL]
input = {Kin, Kout, Imax, IC50, V, Cl}

EQUATION:
; PK
Cc = pkmodel(V, Cl)

; PD - IRM1 (inhibition of production)
ddt_R = Kin * (1 - Imax * Cc / (IC50 + Cc)) - Kout * R

OUTPUT:
output = {Cc, R}
', file.path(MODELS_DIR, "irm1.txt"))

    cat("Monolix model files created in", MODELS_DIR, "\n")
}

# ==============================================================================
# PK Model Benchmarks
# ==============================================================================

benchmark_pk_models <- function() {
    cat("\n======================================================================\n")
    cat("BENCHMARK CATEGORY: PK MODELS\n")
    cat("======================================================================\n")

    # Generate test data
    data <- data.frame(
        ID = rep(1, 49),
        TIME = seq(0, 24, by = 0.5),
        DV = NA,
        AMT = c(100, rep(0, 48)),
        EVID = c(1, rep(0, 48))
    )

    # One-Compartment IV simulation
    tryCatch({
        results <<- rbind(results, run_benchmark(
            function() {
                newProject(modelFile = file.path(MODELS_DIR, "onecomp_iv.txt"))
                setData(dataFile = data)
                setPopulationParameterInformation(V = list(initialValue = 50),
                                                  Cl = list(initialValue = 5))
                runSimulation()
            },
            name = "simulate",
            category = "pk",
            subcategory = "compartmental",
            model = "OneCompIVBolus"
        ))
    }, error = function(e) {
        cat("    SKIPPED:", conditionMessage(e), "\n")
    })

    # Two-Compartment IV simulation
    tryCatch({
        results <<- rbind(results, run_benchmark(
            function() {
                newProject(modelFile = file.path(MODELS_DIR, "twocomp_iv.txt"))
                setData(dataFile = data)
                setPopulationParameterInformation(V1 = list(initialValue = 50),
                                                  V2 = list(initialValue = 100),
                                                  Q = list(initialValue = 10),
                                                  Cl = list(initialValue = 5))
                runSimulation()
            },
            name = "simulate",
            category = "pk",
            subcategory = "compartmental",
            model = "TwoCompIVBolus"
        ))
    }, error = function(e) {
        cat("    SKIPPED:", conditionMessage(e), "\n")
    })

    return(results)
}

# ==============================================================================
# PD Model Benchmarks
# ==============================================================================

benchmark_pd_models <- function() {
    cat("\n======================================================================\n")
    cat("BENCHMARK CATEGORY: PD MODELS\n")
    cat("======================================================================\n")

    # Direct Emax
    tryCatch({
        results <<- rbind(results, run_benchmark(
            function() {
                newProject(modelFile = file.path(MODELS_DIR, "emax_pd.txt"))
                setPopulationParameterInformation(E0 = list(initialValue = 0),
                                                  Emax = list(initialValue = 100),
                                                  EC50 = list(initialValue = 10))
                runSimulation()
            },
            name = "direct_effect",
            category = "pd",
            subcategory = "direct",
            model = "DirectEmax"
        ))
    }, error = function(e) {
        cat("    SKIPPED:", conditionMessage(e), "\n")
    })

    return(results)
}

# ==============================================================================
# Population Benchmarks
# ==============================================================================

benchmark_population <- function() {
    cat("\n======================================================================\n")
    cat("BENCHMARK CATEGORY: POPULATION SIMULATION\n")
    cat("======================================================================\n")

    for (n_subj in c(10, 50, 100, 500, 1000)) {
        tryCatch({
            # Generate population data
            data <- data.frame(
                ID = rep(1:n_subj, each = 49),
                TIME = rep(seq(0, 24, by = 0.5), n_subj),
                DV = NA,
                AMT = rep(c(100, rep(0, 48)), n_subj),
                EVID = rep(c(1, rep(0, 48)), n_subj)
            )

            results <<- rbind(results, run_benchmark(
                function() {
                    newProject(modelFile = file.path(MODELS_DIR, "twocomp_iv.txt"))
                    setData(dataFile = data)
                    setPopulationParameterInformation(
                        V1 = list(initialValue = 50),
                        V2 = list(initialValue = 100),
                        Q = list(initialValue = 10),
                        Cl = list(initialValue = 5)
                    )
                    # Add IIV
                    setIndividualParameterVariability(
                        V1 = TRUE, V2 = TRUE, Q = TRUE, Cl = TRUE
                    )
                    runSimulation()
                },
                name = "population_simulation",
                category = "population",
                subcategory = "iiv",
                model = "TwoCompIV_IIV",
                n_subjects = n_subj
            ))
        }, error = function(e) {
            cat(sprintf("    SKIPPED (n=%d): %s\n", n_subj, conditionMessage(e)))
        })
    }

    return(results)
}

# ==============================================================================
# Main Execution
# ==============================================================================

cat("Creating Monolix model files...\n")
create_monolix_models()

results <- benchmark_pk_models()
results <- benchmark_pd_models()
results <- benchmark_population()

# Save results
output_file <- file.path(RESULTS_DIR, "monolix_comprehensive_benchmarks.csv")
write.csv(results, output_file, row.names = FALSE)

cat("\n==============================================================================\n")
cat("MONOLIX BENCHMARK COMPLETE\n")
cat(sprintf("Results saved to: %s\n", output_file))
cat(sprintf("Total benchmarks: %d\n", nrow(results)))
cat("==============================================================================\n")
