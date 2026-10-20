#!/usr/bin/env Rscript
#' NeoPKPD Benchmark Suite - mrgsolve Comparison
#' ==================================================
#'
#' Benchmarks mrgsolve for comparison with NeoPKPD.
#' mrgsolve is a popular R package for ODE-based PK/PD simulation.
#'
#' Installation:
#'   install.packages("mrgsolve")
#'   install.packages("dplyr")
#'   install.packages("microbenchmark")
#'
#' Usage:
#'   Rscript benchmark_mrgsolve.R
#'
#' Output:
#'   benchmarks/results/mrgsolve_benchmarks.csv

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
  library(dplyr)
  library(microbenchmark)
})

# Configuration
# Get script directory robustly
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("--file=", "", file_arg))))
  }
  return(getwd())
}

script_dir <- get_script_dir()
CONFIG <- list(
  n_runs = 100,
  n_warmup = 5,
  random_seed = 12345,
  output_dir = file.path(dirname(script_dir), "results")
)

# Create output directory
dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat(strrep("=", 70), "\n")
cat("mrgsolve Benchmark Suite\n")
cat(strrep("=", 70), "\n")
cat("R Version:", R.version.string, "\n")
cat("mrgsolve Version:", as.character(packageVersion("mrgsolve")), "\n")
cat("Runs per benchmark:", CONFIG$n_runs, "\n")
cat("Timestamp:", format(Sys.time()), "\n")
cat(strrep("=", 70), "\n\n")

set.seed(CONFIG$random_seed)

# =============================================================================
# Model Definitions (mrgsolve format)
# =============================================================================

# One-compartment IV bolus
mod_1comp <- mcode("onecomp_iv", '
$PARAM CL = 5, V = 50
$CMT CENT
$ODE
dxdt_CENT = -CL/V * CENT;
$TABLE
double CP = CENT/V;
$CAPTURE CP
', compile = TRUE)

# Two-compartment IV bolus
mod_2comp <- mcode("twocomp_iv", '
$PARAM CL = 5, V1 = 50, V2 = 100, Q = 10
$CMT CENT PERIPH
$ODE
double C1 = CENT/V1;
double C2 = PERIPH/V2;
dxdt_CENT = -CL*C1 - Q*C1 + Q*C2;
dxdt_PERIPH = Q*C1 - Q*C2;
$TABLE
double CP = CENT/V1;
$CAPTURE CP
', compile = TRUE)

# Two-compartment oral
mod_2comp_oral <- mcode("twocomp_oral", '
$PARAM CL = 5, V1 = 50, V2 = 100, Q = 10, KA = 1.5
$CMT GUT CENT PERIPH
$ODE
double C1 = CENT/V1;
double C2 = PERIPH/V2;
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - CL*C1 - Q*C1 + Q*C2;
dxdt_PERIPH = Q*C1 - Q*C2;
$TABLE
double CP = CENT/V1;
$CAPTURE CP
', compile = TRUE)

# Michaelis-Menten elimination
mod_mm <- mcode("mm_elim", '
$PARAM VMAX = 50, KM = 10, V = 50
$CMT CENT
$ODE
double C = CENT/V;
dxdt_CENT = -VMAX * C / (KM + C);
$TABLE
double CP = CENT/V;
$CAPTURE CP
', compile = TRUE)

# Direct Emax PKPD
mod_emax <- mcode("pkpd_emax", '
$PARAM CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10
$CMT CENT
$ODE
dxdt_CENT = -CL/V * CENT;
$TABLE
double CP = CENT/V;
double EFFECT = E0 + EMAX * CP / (EC50 + CP);
$CAPTURE CP EFFECT
', compile = TRUE)

# Indirect response model (IRM3 - inhibition of kout)
mod_irm <- mcode("pkpd_irm", '
$PARAM CL = 5, V = 50, R0 = 100, KOUT = 0.1, IMAX = 0.9, IC50 = 10
$CMT CENT RESP
$MAIN
RESP_0 = R0;
$ODE
double CP = CENT/V;
double KIN = R0 * KOUT;
double INH = 1 - IMAX * CP / (IC50 + CP);
dxdt_CENT = -CL/V * CENT;
dxdt_RESP = KIN - KOUT * INH * RESP;
$TABLE
double EFFECT = RESP;
$CAPTURE CP EFFECT
', compile = TRUE)

# =============================================================================
# Benchmark Functions
# =============================================================================

run_benchmark <- function(name, category, model_name, n_subjects, expr) {
  cat("  Running:", name, "(", model_name, ")...\n")

  # Warmup
  for (i in 1:CONFIG$n_warmup) {
    eval(expr)
  }

  # Timed runs
  mb <- microbenchmark(
    eval(expr),
    times = CONFIG$n_runs,
    unit = "ms"
  )

  times_ms <- mb$time / 1e6  # nanoseconds to milliseconds

  result <- data.frame(
    category = category,
    name = name,
    model = model_name,
    n_subjects = n_subjects,
    n_runs = CONFIG$n_runs,
    mean_ms = mean(times_ms),
    std_ms = sd(times_ms),
    median_ms = median(times_ms),
    min_ms = min(times_ms),
    max_ms = max(times_ms),
    ci_lower_ms = quantile(times_ms, 0.025),
    ci_upper_ms = quantile(times_ms, 0.975),
    timestamp = format(Sys.time()),
    r_version = R.version.string,
    mrgsolve_version = as.character(packageVersion("mrgsolve")),
    stringsAsFactors = FALSE
  )

  cat(sprintf("    Mean: %.3f ms (±%.3f), Median: %.3f ms\n",
              result$mean_ms, result$std_ms, result$median_ms))

  return(result)
}

# =============================================================================
# Benchmark Execution
# =============================================================================

results <- data.frame()

# Simulation parameters
dose_data <- data.frame(ID = 1, time = 0, amt = 100, cmt = 1, evid = 1)
sim_times <- seq(0, 24, by = 0.5)

# -----------------------------------------------------------------------------
# 1. Single Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: Single Simulation Speed\n")
cat(strrep("=", 60), "\n")

# 1a. One-compartment IV
results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "OneCompIV", 1,
  quote({
    mod_1comp %>%
      data_set(dose_data) %>%
      mrgsim(end = 24, delta = 0.5)
  })
))

# 1b. Two-compartment IV
results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "TwoCompIV", 1,
  quote({
    mod_2comp %>%
      data_set(dose_data) %>%
      mrgsim(end = 24, delta = 0.5)
  })
))

# 1c. Two-compartment Oral
dose_oral <- data.frame(ID = 1, time = 0, amt = 100, cmt = 1, evid = 1)
results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "TwoCompOral", 1,
  quote({
    mod_2comp_oral %>%
      data_set(dose_oral) %>%
      mrgsim(end = 24, delta = 0.5)
  })
))

# 1d. Michaelis-Menten
results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "MichaelisMenten", 1,
  quote({
    mod_mm %>%
      data_set(dose_data) %>%
      mrgsim(end = 24, delta = 0.5)
  })
))

# -----------------------------------------------------------------------------
# 2. Population Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: Population Simulation Speed\n")
cat(strrep("=", 60), "\n")

# Create population models with IIV
mod_2comp_pop <- mcode("twocomp_pop", '
$PARAM CL = 5, V1 = 50, V2 = 100, Q = 10
$OMEGA 0.09 0.04  // 30% CV CL, 20% CV V1
$CMT CENT PERIPH
$MAIN
double CLi = CL * exp(ETA(1));
double V1i = V1 * exp(ETA(2));
$ODE
double C1 = CENT/V1i;
double C2 = PERIPH/V2;
dxdt_CENT = -CLi*C1 - Q*C1 + Q*C2;
dxdt_PERIPH = Q*C1 - Q*C2;
$TABLE
double CP = CENT/V1i;
$CAPTURE CP
', compile = TRUE)

for (n_subj in c(50, 100, 500, 1000)) {
  pop_dose <- expand.grid(ID = 1:n_subj) %>%
    mutate(time = 0, amt = 100, cmt = 1, evid = 1)

  results <- rbind(results, run_benchmark(
    "population_simulation", "population", "TwoCompIV", n_subj,
    bquote({
      mod_2comp_pop %>%
        data_set(.(pop_dose)) %>%
        mrgsim(end = 24, delta = 0.5)
    })
  ))
}

# -----------------------------------------------------------------------------
# 3. PKPD Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: PKPD Coupled Simulation Speed\n")
cat(strrep("=", 60), "\n")

# 3a. Direct Emax
results <- rbind(results, run_benchmark(
  "pkpd_simulation", "pkpd", "DirectEmax", 1,
  quote({
    mod_emax %>%
      data_set(dose_data) %>%
      mrgsim(end = 72, delta = 1)
  })
))

# 3b. Indirect Response
results <- rbind(results, run_benchmark(
  "pkpd_simulation", "pkpd", "IndirectResponse", 1,
  quote({
    mod_irm %>%
      data_set(dose_data) %>%
      mrgsim(end = 72, delta = 1)
  })
))

# =============================================================================
# Save Results
# =============================================================================

output_file <- file.path(CONFIG$output_dir, "mrgsolve_benchmarks.csv")
write.csv(results, output_file, row.names = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK COMPLETE\n")
cat("Results saved to:", output_file, "\n")
cat(strrep("=", 70), "\n")
