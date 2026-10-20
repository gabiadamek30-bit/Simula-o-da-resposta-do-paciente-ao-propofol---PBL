#!/usr/bin/env Rscript
#' NeoPKPD Benchmark Suite - nlmixr2 Comparison
#' ==================================================
#'
#' Benchmarks nlmixr2 for comparison with NeoPKPD.
#' nlmixr2 is an R package for nonlinear mixed-effects modeling.
#'
#' Installation:
#'   install.packages("nlmixr2")
#'   install.packages("rxode2")
#'   install.packages("dplyr")
#'   install.packages("microbenchmark")
#'
#' Usage:
#'   Rscript benchmark_nlmixr2.R
#'
#' Output:
#'   benchmarks/results/nlmixr2_benchmarks.csv

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
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
cat("nlmixr2 Benchmark Suite\n")
cat(strrep("=", 70), "\n")
cat("R Version:", R.version.string, "\n")
cat("nlmixr2 Version:", as.character(packageVersion("nlmixr2")), "\n")
cat("rxode2 Version:", as.character(packageVersion("rxode2")), "\n")
cat("Runs per benchmark:", CONFIG$n_runs, "\n")
cat("Timestamp:", format(Sys.time()), "\n")
cat(strrep("=", 70), "\n\n")

set.seed(CONFIG$random_seed)

# =============================================================================
# Model Definitions (rxode2 format)
# =============================================================================

# One-compartment IV bolus
mod_1comp <- rxode2({
  C <- centr / V
  d/dt(centr) <- -CL/V * centr
})

# Two-compartment IV bolus
mod_2comp <- rxode2({
  C1 <- centr / V1
  C2 <- periph / V2
  d/dt(centr) <- -CL * C1 - Q * C1 + Q * C2
  d/dt(periph) <- Q * C1 - Q * C2
  C <- C1
})

# Two-compartment oral
mod_2comp_oral <- rxode2({
  C1 <- centr / V1
  C2 <- periph / V2
  d/dt(depot) <- -KA * depot
  d/dt(centr) <- KA * depot - CL * C1 - Q * C1 + Q * C2
  d/dt(periph) <- Q * C1 - Q * C2
  C <- C1
})

# Michaelis-Menten elimination
mod_mm <- rxode2({
  C <- centr / V
  d/dt(centr) <- -VMAX * C / (KM + C)
})

# Direct Emax PKPD
mod_emax <- rxode2({
  C <- centr / V
  d/dt(centr) <- -CL/V * centr
  EFFECT <- E0 + EMAX * C / (EC50 + C)
})

# Indirect response model (IRM3)
mod_irm <- rxode2({
  C <- centr / V
  KIN <- R0 * KOUT
  INH <- 1 - IMAX * C / (IC50 + C)
  d/dt(centr) <- -CL/V * centr
  d/dt(resp) <- KIN - KOUT * INH * resp
  EFFECT <- resp
})

# =============================================================================
# Benchmark Functions
# =============================================================================

run_benchmark <- function(name, category, model_name, n_subjects, expr) {
  cat("  Running:", name, "(", model_name, ")...\n")

  # Warmup
  for (i in 1:CONFIG$n_warmup) {
    tryCatch(eval(expr), error = function(e) NULL)
  }

  # Timed runs
  times_ms <- numeric(CONFIG$n_runs)
  for (i in 1:CONFIG$n_runs) {
    start_time <- Sys.time()
    tryCatch(eval(expr), error = function(e) NULL)
    end_time <- Sys.time()
    times_ms[i] <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
  }

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
    nlmixr2_version = as.character(packageVersion("nlmixr2")),
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
sim_times <- seq(0, 24, by = 0.5)

# -----------------------------------------------------------------------------
# 1. Single Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: Single Simulation Speed\n")
cat(strrep("=", 60), "\n")

# Event table for dosing
ev_1comp <- et(amt = 100, cmt = 1) %>% et(seq(0, 24, by = 0.5))
params_1comp <- c(CL = 5, V = 50)

# 1a. One-compartment IV
results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "OneCompIV", 1,
  quote({
    rxSolve(mod_1comp, params_1comp, ev_1comp)
  })
))

# 1b. Two-compartment IV
ev_2comp <- et(amt = 100, cmt = 1) %>% et(seq(0, 24, by = 0.5))
params_2comp <- c(CL = 5, V1 = 50, V2 = 100, Q = 10)

results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "TwoCompIV", 1,
  quote({
    rxSolve(mod_2comp, params_2comp, ev_2comp)
  })
))

# 1c. Two-compartment Oral
params_2comp_oral <- c(CL = 5, V1 = 50, V2 = 100, Q = 10, KA = 1.5)

results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "TwoCompOral", 1,
  quote({
    rxSolve(mod_2comp_oral, params_2comp_oral, ev_2comp)
  })
))

# 1d. Michaelis-Menten
params_mm <- c(VMAX = 50, KM = 10, V = 50)

results <- rbind(results, run_benchmark(
  "single_simulation", "simulation", "MichaelisMenten", 1,
  quote({
    rxSolve(mod_mm, params_mm, ev_1comp)
  })
))

# -----------------------------------------------------------------------------
# 2. Population Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: Population Simulation Speed\n")
cat(strrep("=", 60), "\n")

# Population model with IIV
mod_2comp_pop <- rxode2({
  CLi <- CL * exp(eta.CL)
  V1i <- V1 * exp(eta.V1)
  C1 <- centr / V1i
  C2 <- periph / V2
  d/dt(centr) <- -CLi * C1 - Q * C1 + Q * C2
  d/dt(periph) <- Q * C1 - Q * C2
  C <- C1
})

omega <- lotri({
  eta.CL ~ 0.09  # 30% CV
  eta.V1 ~ 0.04  # 20% CV
})

for (n_subj in c(50, 100, 500, 1000)) {
  # Create event table for population
  ev_pop <- et(amt = 100, cmt = 1) %>%
    et(seq(0, 24, by = 0.5)) %>%
    et(id = 1:n_subj)

  results <- rbind(results, run_benchmark(
    "population_simulation", "population", "TwoCompIV", n_subj,
    bquote({
      rxSolve(mod_2comp_pop, params_2comp, .(ev_pop),
              omega = omega, nSub = .(n_subj))
    })
  ))
}

# -----------------------------------------------------------------------------
# 3. PKPD Simulation Benchmarks
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: PKPD Coupled Simulation Speed\n")
cat(strrep("=", 60), "\n")

ev_pkpd <- et(amt = 100, cmt = 1) %>% et(seq(0, 72, by = 1))

# 3a. Direct Emax
params_emax <- c(CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10)

results <- rbind(results, run_benchmark(
  "pkpd_simulation", "pkpd", "DirectEmax", 1,
  quote({
    rxSolve(mod_emax, params_emax, ev_pkpd)
  })
))

# 3b. Indirect Response
params_irm <- c(CL = 5, V = 50, R0 = 100, KOUT = 0.1, IMAX = 0.9, IC50 = 10)
ev_irm <- et(amt = 100, cmt = 1) %>% et(seq(0, 72, by = 1))
# Set initial response
inits_irm <- c(resp = 100)

results <- rbind(results, run_benchmark(
  "pkpd_simulation", "pkpd", "IndirectResponse", 1,
  quote({
    rxSolve(mod_irm, params_irm, ev_irm, inits = inits_irm)
  })
))

# -----------------------------------------------------------------------------
# 4. Estimation Benchmarks (if data available)
# -----------------------------------------------------------------------------

cat("\n", strrep("=", 60), "\n")
cat("BENCHMARK: Parameter Estimation Speed\n")
cat(strrep("=", 60), "\n")

# Generate synthetic data for estimation
generate_data <- function(n_subj, n_obs_per_subj) {
  set.seed(CONFIG$random_seed)

  data <- data.frame()
  for (id in 1:n_subj) {
    # Individual parameters with IIV
    eta_cl <- rnorm(1, 0, 0.3)
    eta_v <- rnorm(1, 0, 0.2)
    cli <- 5 * exp(eta_cl)
    vi <- 50 * exp(eta_v)

    # Observation times
    times <- c(0.5, 1, 2, 4, 8, 12, 24)

    # Simulate concentrations (one-comp IV bolus)
    conc <- (100/vi) * exp(-cli/vi * times)

    # Add residual error (proportional)
    conc <- conc * exp(rnorm(length(times), 0, 0.1))

    subj_data <- data.frame(
      ID = id,
      TIME = c(0, times),
      DV = c(NA, conc),
      AMT = c(100, rep(0, length(times))),
      EVID = c(1, rep(0, length(times))),
      CMT = 1
    )
    data <- rbind(data, subj_data)
  }
  return(data)
}

# nlmixr2 model for estimation
one.cmt <- function() {
  ini({
    tCL <- log(5)
    tV <- log(50)
    eta.CL ~ 0.09
    eta.V ~ 0.04
    prop.err <- 0.1
  })
  model({
    CL <- exp(tCL + eta.CL)
    V <- exp(tV + eta.V)
    d/dt(centr) <- -CL/V * centr
    cp <- centr/V
    cp ~ prop(prop.err)
  })
}

# Only run estimation benchmarks with smaller dataset (time-consuming)
n_subj_est <- 50
est_data <- generate_data(n_subj_est, 7)

cat("  Note: Estimation benchmarks use N=", n_subj_est, " subjects\n")
cat("  Running FOCE-I estimation (this may take several minutes)...\n")

# Reduced runs for estimation (it's slow)
est_times <- numeric(10)  # Only 10 runs for estimation
for (i in 1:10) {
  start_time <- Sys.time()
  tryCatch({
    fit <- nlmixr2(one.cmt, est_data, est = "focei",
                   control = foceiControl(print = 0))
  }, error = function(e) NULL)
  end_time <- Sys.time()
  est_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs")) * 1000
}

results <- rbind(results, data.frame(
  category = "estimation",
  name = "estimation_focei",
  model = "OneCompIV",
  n_subjects = n_subj_est,
  n_runs = 10,
  mean_ms = mean(est_times),
  std_ms = sd(est_times),
  median_ms = median(est_times),
  min_ms = min(est_times),
  max_ms = max(est_times),
  ci_lower_ms = quantile(est_times, 0.025),
  ci_upper_ms = quantile(est_times, 0.975),
  timestamp = format(Sys.time()),
  r_version = R.version.string,
  nlmixr2_version = as.character(packageVersion("nlmixr2")),
  stringsAsFactors = FALSE
))

cat(sprintf("    Mean: %.1f ms (±%.1f), Median: %.1f ms\n",
            mean(est_times), sd(est_times), median(est_times)))

# =============================================================================
# Save Results
# =============================================================================

output_file <- file.path(CONFIG$output_dir, "nlmixr2_benchmarks.csv")
write.csv(results, output_file, row.names = FALSE)

cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK COMPLETE\n")
cat("Results saved to:", output_file, "\n")
cat(strrep("=", 70), "\n")
