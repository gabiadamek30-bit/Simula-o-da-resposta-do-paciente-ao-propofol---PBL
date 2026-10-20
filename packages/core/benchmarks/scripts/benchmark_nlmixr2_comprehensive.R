#!/usr/bin/env Rscript
#' NeoPKPD Benchmark Suite - nlmixr2 Comprehensive Comparison
#' ==============================================================
#'
#' Comprehensive nlmixr2/rxode2 benchmarks matching NeoPKPD categories.
#'
#' Usage:
#'   Rscript benchmark_nlmixr2_comprehensive.R

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(nlmixr2)
  library(rxode2)
  library(dplyr)
  library(microbenchmark)
})

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
  n_runs = 50,
  n_warmup = 3,
  random_seed = 12345,
  output_dir = file.path(dirname(script_dir), "results")
)

dir.create(CONFIG$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat(strrep("=", 80), "\n")
cat("nlmixr2/rxode2 COMPREHENSIVE Benchmark Suite\n")
cat(strrep("=", 80), "\n")
cat("R Version:", R.version.string, "\n")
cat("nlmixr2 Version:", as.character(packageVersion("nlmixr2")), "\n")
cat("rxode2 Version:", as.character(packageVersion("rxode2")), "\n")
cat("Runs per benchmark:", CONFIG$n_runs, "\n")
cat("Timestamp:", format(Sys.time()), "\n")
cat(strrep("=", 80), "\n\n")

set.seed(CONFIG$random_seed)

# =============================================================================
# Benchmark Infrastructure
# =============================================================================

run_benchmark <- function(expr, name, category, subcategory = "", model, n_subjects = 1) {
  subcat_str <- ifelse(subcategory == "", "", paste0(" [", subcategory, "]"))
  cat("  Running:", name, "(", model, ")", subcat_str, "...\n")

  # Warmup
  for (i in 1:CONFIG$n_warmup) {
    tryCatch(eval(expr), error = function(e) NULL)
  }

  # Timed runs
  mb <- tryCatch({
    microbenchmark(eval(expr), times = CONFIG$n_runs, unit = "ms")
  }, error = function(e) {
    cat("    ERROR:", e$message, "\n")
    return(NULL)
  })

  if (is.null(mb)) return(NULL)

  times_ms <- mb$time / 1e6

  cat(sprintf("    Mean: %.3f ms (±%.3f), Median: %.3f ms\n",
              mean(times_ms), sd(times_ms), median(times_ms)))

  data.frame(
    category = category,
    subcategory = subcategory,
    name = name,
    model = model,
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
    rxode2_version = as.character(packageVersion("rxode2"))
  )
}

# =============================================================================
# Model Definitions
# =============================================================================

cat("Building rxode2 models...\n")

# 1-compartment IV
onecomp_iv <- rxode2({
  d/dt(CENT) = -CL/V * CENT
  CP = CENT/V
})

# 1-compartment oral
onecomp_oral <- rxode2({
  d/dt(GUT) = -KA * GUT
  d/dt(CENT) = KA * GUT - CL/V * CENT
  CP = CENT/V
})

# 2-compartment IV
twocomp_iv <- rxode2({
  k10 = CL/V1
  k12 = Q/V1
  k21 = Q/V2
  d/dt(CENT) = -k10*CENT - k12*CENT + k21*PERIPH
  d/dt(PERIPH) = k12*CENT - k21*PERIPH
  CP = CENT/V1
})

# 2-compartment oral
twocomp_oral <- rxode2({
  k10 = CL/V1
  k12 = Q/V1
  k21 = Q/V2
  d/dt(GUT) = -KA * GUT
  d/dt(CENT) = KA*GUT - k10*CENT - k12*CENT + k21*PERIPH
  d/dt(PERIPH) = k12*CENT - k21*PERIPH
  CP = CENT/V1
})

# 3-compartment IV
threecomp_iv <- rxode2({
  k10 = CL/V1
  k12 = Q2/V1
  k21 = Q2/V2
  k13 = Q3/V1
  k31 = Q3/V3
  d/dt(CENT) = -k10*CENT - k12*CENT + k21*PERIPH1 - k13*CENT + k31*PERIPH2
  d/dt(PERIPH1) = k12*CENT - k21*PERIPH1
  d/dt(PERIPH2) = k13*CENT - k31*PERIPH2
  CP = CENT/V1
})

# Michaelis-Menten elimination
mm_elim <- rxode2({
  C = CENT/V
  d/dt(CENT) = -VMAX * C / (KM + C)
  CP = CENT/V
})

# Direct Emax PD
pkpd_emax <- rxode2({
  d/dt(CENT) = -CL/V * CENT
  CP = CENT/V
  EFF = E0 + EMAX * CP / (EC50 + CP)
})

# Sigmoid Emax PD
pkpd_sigmoid <- rxode2({
  d/dt(CENT) = -CL/V * CENT
  CP = CENT/V
  EFF = E0 + EMAX * CP^GAMMA / (EC50^GAMMA + CP^GAMMA)
})

# IRM1
pkpd_irm1 <- rxode2({
  d/dt(CENT) = -CL/V * CENT
  CP = CENT/V
  INH = IMAX * CP / (IC50 + CP)
  d/dt(RESP) = KIN * (1 - INH) - KOUT * RESP
})

# Biophase equilibration
pkpd_biophase <- rxode2({
  d/dt(CENT) = -CL/V * CENT
  CP = CENT/V
  d/dt(CE) = KE0 * (CP - CE)
  EFF = E0 + EMAX * CE / (EC50 + CE)
})

cat("All models compiled successfully.\n\n")

# =============================================================================
# Benchmark Execution
# =============================================================================

all_results <- data.frame()

# Common parameters
params_1comp <- c(CL = 5, V = 50)
params_1comp_oral <- c(KA = 1.5, CL = 5, V = 50)
params_2comp <- c(CL = 5, V1 = 50, Q = 10, V2 = 100)
params_2comp_oral <- c(KA = 1.5, CL = 5, V1 = 50, Q = 10, V2 = 100)
params_3comp <- c(CL = 5, V1 = 50, Q2 = 8, V2 = 80, Q3 = 3, V3 = 150)
params_mm <- c(VMAX = 50, KM = 10, V = 50)
params_emax <- c(CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10)
params_sigmoid <- c(CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10, GAMMA = 2)
params_irm1 <- c(CL = 5, V = 50, KIN = 10, KOUT = 0.1, IMAX = 0.9, IC50 = 5)
params_biophase <- c(CL = 5, V = 50, KE0 = 0.5, E0 = 0, EMAX = 100, EC50 = 10)

# Event table
ev <- et(amt = 100, cmt = 1) %>%
  et(seq(0, 24, by = 0.5))

ev_72h <- et(amt = 100, cmt = 1) %>%
  et(seq(0, 72, by = 1))

# Multiple dosing
ev_multi <- et(amt = 100, cmt = 1, ii = 24, addl = 6) %>%
  et(seq(0, 168, by = 2))

# -----------------------------------------------------------------------------
# 1. PK Models
# -----------------------------------------------------------------------------
cat(strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PK MODELS\n")
cat(strrep("=", 70), "\n")

# One-compartment IV
result <- run_benchmark(
  quote(rxSolve(onecomp_iv, params_1comp, ev)),
  "simulate", "pk", "compartmental", "OneCompIVBolus"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# One-compartment oral
result <- run_benchmark(
  quote(rxSolve(onecomp_oral, params_1comp_oral, ev)),
  "simulate", "pk", "compartmental", "OneCompOralFirstOrder"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Two-compartment IV
result <- run_benchmark(
  quote(rxSolve(twocomp_iv, params_2comp, ev)),
  "simulate", "pk", "compartmental", "TwoCompIVBolus"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Two-compartment oral
result <- run_benchmark(
  quote(rxSolve(twocomp_oral, params_2comp_oral, ev)),
  "simulate", "pk", "compartmental", "TwoCompOral"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Three-compartment IV
result <- run_benchmark(
  quote(rxSolve(threecomp_iv, params_3comp, ev)),
  "simulate", "pk", "compartmental", "ThreeCompIVBolus"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Michaelis-Menten
result <- run_benchmark(
  quote(rxSolve(mm_elim, params_mm, ev)),
  "simulate", "pk", "nonlinear", "MichaelisMenten"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Multiple dosing (7 days)
result <- run_benchmark(
  quote(rxSolve(twocomp_iv, params_2comp, ev_multi)),
  "simulate_multiple_dose", "pk", "dosing", "TwoCompIV_7days"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# -----------------------------------------------------------------------------
# 2. PD Models
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PD MODELS\n")
cat(strrep("=", 70), "\n")

# Direct Emax
result <- run_benchmark(
  quote(rxSolve(pkpd_emax, params_emax, ev)),
  "direct_effect", "pd", "direct", "DirectEmax"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Sigmoid Emax
result <- run_benchmark(
  quote(rxSolve(pkpd_sigmoid, params_sigmoid, ev)),
  "direct_effect", "pd", "direct", "SigmoidEmax"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# IRM1 (need initial condition for RESP)
params_irm1_full <- c(params_irm1, RESP = 100)  # Baseline = Kin/Kout
result <- run_benchmark(
  quote(rxSolve(pkpd_irm1, params_irm1_full, ev_72h)),
  "indirect_response", "pd", "indirect", "IRM1"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Biophase
result <- run_benchmark(
  quote(rxSolve(pkpd_biophase, params_biophase, ev_72h)),
  "effect_compartment", "pd", "biophase", "BiophaseEquilibration"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# -----------------------------------------------------------------------------
# 3. PKPD Coupled
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PKPD COUPLED\n")
cat(strrep("=", 70), "\n")

result <- run_benchmark(
  quote(rxSolve(pkpd_emax, params_emax, ev_72h)),
  "pkpd_coupled", "pkpd", "direct", "OneComp_DirectEmax"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

result <- run_benchmark(
  quote(rxSolve(pkpd_sigmoid, params_sigmoid, ev_72h)),
  "pkpd_coupled", "pkpd", "direct", "OneComp_SigmoidEmax"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# -----------------------------------------------------------------------------
# 4. NCA
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: NCA\n")
cat(strrep("=", 70), "\n")

# Generate PK data
sim_result <- rxSolve(twocomp_iv, params_2comp, ev)
pk_data <- as.data.frame(sim_result)
times_nca <- pk_data$time
conc_nca <- pk_data$CP

# Cmax/Tmax
result <- run_benchmark(
  quote({
    cmax <- max(conc_nca)
    tmax <- times_nca[which.max(conc_nca)]
    c(cmax, tmax)
  }),
  "cmax_tmax", "nca", "exposure", "Cmax_Tmax"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# AUC
result <- run_benchmark(
  quote({
    n <- length(times_nca)
    auc <- sum((conc_nca[-1] + conc_nca[-n]) / 2 * diff(times_nca))
    auc
  }),
  "auc_0_t", "nca", "exposure", "AUC_0_t"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Lambda-z
result <- run_benchmark(
  quote({
    n <- length(conc_nca)
    idx <- (n-4):n
    lm_fit <- lm(log(conc_nca[idx] + 1e-10) ~ times_nca[idx])
    lambda_z <- -coef(lm_fit)[2]
    lambda_z
  }),
  "lambda_z", "nca", "kinetics", "Lambda_z"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# -----------------------------------------------------------------------------
# 5. Population Simulation
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: POPULATION SIMULATION\n")
cat(strrep("=", 70), "\n")

# Population sizes
pop_sizes <- c(10, 50, 100, 500, 1000)

for (n_subj in pop_sizes) {
  # Generate population with IIV
  omega <- matrix(c(0.09, 0, 0, 0.04), 2, 2)
  eta <- MASS::mvrnorm(n_subj, c(0, 0), omega)

  # Create individual parameters
  pop_params <- data.frame(
    id = 1:n_subj,
    CL = 5 * exp(eta[, 1]),
    V1 = 50 * exp(eta[, 2]),
    Q = 10,
    V2 = 100
  )

  ev_pop <- et(amt = 100, cmt = 1) %>%
    et(seq(0, 24, by = 0.5))

  result <- run_benchmark(
    quote(rxSolve(twocomp_iv, pop_params, ev_pop)),
    "population_simulation", "population", "iiv",
    "TwoCompIV_IIV", n_subjects = n_subj
  )
  if (!is.null(result)) all_results <- rbind(all_results, result)
}

# -----------------------------------------------------------------------------
# 6. Sensitivity Analysis
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: SENSITIVITY ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Local sensitivity
result <- run_benchmark(
  quote({
    base <- rxSolve(onecomp_iv, params_1comp, ev)

    params_pert <- params_1comp
    params_pert["CL"] <- 5.5
    pert <- rxSolve(onecomp_iv, params_pert, ev)

    sens <- (as.data.frame(pert)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP
    sens
  }),
  "local_single", "sensitivity", "local", "LocalSingle"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# Multi-parameter sensitivity
result <- run_benchmark(
  quote({
    base <- rxSolve(onecomp_iv, params_1comp, ev)

    params_cl <- params_1comp
    params_cl["CL"] <- 5.5
    pert_cl <- rxSolve(onecomp_iv, params_cl, ev)

    params_v <- params_1comp
    params_v["V"] <- 55
    pert_v <- rxSolve(onecomp_iv, params_v, ev)

    list(
      sens_cl = (as.data.frame(pert_cl)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP,
      sens_v = (as.data.frame(pert_v)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP
    )
  }),
  "local_multi", "sensitivity", "local", "LocalMulti"
)
if (!is.null(result)) all_results <- rbind(all_results, result)

# =============================================================================
# Save Results
# =============================================================================

if (nrow(all_results) > 0) {
  output_file <- file.path(CONFIG$output_dir, "nlmixr2_comprehensive_benchmarks.csv")
  write.csv(all_results, output_file, row.names = FALSE)

  cat("\n", strrep("=", 80), "\n")
  cat("COMPREHENSIVE BENCHMARK COMPLETE\n")
  cat("Results saved to:", output_file, "\n")
  cat("Total benchmarks:", nrow(all_results), "\n")
  cat(strrep("=", 80), "\n")

  # Print summary
  cat("\n", strrep("-", 80), "\n")
  cat("SUMMARY BY CATEGORY\n")
  cat(strrep("-", 80), "\n")

  summary_df <- all_results %>%
    group_by(category) %>%
    summarize(
      n_benchmarks = n(),
      avg_ms = mean(mean_ms),
      .groups = "drop"
    )

  for (i in 1:nrow(summary_df)) {
    cat(sprintf("%-20s: %3d benchmarks, avg %.3f ms\n",
                summary_df$category[i],
                summary_df$n_benchmarks[i],
                summary_df$avg_ms[i]))
  }
  cat(strrep("-", 80), "\n")
} else {
  cat("\nNo benchmark results collected.\n")
}
