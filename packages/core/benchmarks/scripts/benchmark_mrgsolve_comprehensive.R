#!/usr/bin/env Rscript
#' NeoPKPD Benchmark Suite - mrgsolve Comprehensive Comparison
#' ==============================================================
#'
#' Comprehensive mrgsolve benchmarks matching NeoPKPD categories:
#' - PK Models (1-comp, 2-comp, 3-comp, oral, MM)
#' - PD Models (Emax, Hill, IRM)
#' - PKPD Coupled
#' - NCA calculations
#' - Population simulation
#' - Sensitivity analysis
#'
#' Usage:
#'   Rscript benchmark_mrgsolve_comprehensive.R

# =============================================================================
# Setup
# =============================================================================

suppressPackageStartupMessages({
  library(mrgsolve)
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
cat("mrgsolve COMPREHENSIVE Benchmark Suite\n")
cat(strrep("=", 80), "\n")
cat("R Version:", R.version.string, "\n")
cat("mrgsolve Version:", as.character(packageVersion("mrgsolve")), "\n")
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
  mb <- microbenchmark(eval(expr), times = CONFIG$n_runs, unit = "ms")
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
    mrgsolve_version = as.character(packageVersion("mrgsolve"))
  )
}

# =============================================================================
# Model Definitions
# =============================================================================

cat("Building models...\n")

# 1-compartment IV
onecomp_iv <- mcode("onecomp_iv", "
$PARAM CL = 5, V = 50
$CMT CENT
$ODE dxdt_CENT = -CL/V * CENT;
$TABLE double CP = CENT/V;
$CAPTURE CP
", compile = TRUE)

# 1-compartment oral
onecomp_oral <- mcode("onecomp_oral", "
$PARAM KA = 1.5, CL = 5, V = 50
$CMT GUT CENT
$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - CL/V * CENT;
$TABLE double CP = CENT/V;
$CAPTURE CP
", compile = TRUE)

# 2-compartment IV
twocomp_iv <- mcode("twocomp_iv", "
$PARAM CL = 5, V1 = 50, Q = 10, V2 = 100
$CMT CENT PERIPH
$ODE
double k10 = CL/V1;
double k12 = Q/V1;
double k21 = Q/V2;
dxdt_CENT = -k10*CENT - k12*CENT + k21*PERIPH;
dxdt_PERIPH = k12*CENT - k21*PERIPH;
$TABLE double CP = CENT/V1;
$CAPTURE CP
", compile = TRUE)

# 2-compartment oral
twocomp_oral <- mcode("twocomp_oral", "
$PARAM KA = 1.5, CL = 5, V1 = 50, Q = 10, V2 = 100
$CMT GUT CENT PERIPH
$ODE
double k10 = CL/V1;
double k12 = Q/V1;
double k21 = Q/V2;
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA*GUT - k10*CENT - k12*CENT + k21*PERIPH;
dxdt_PERIPH = k12*CENT - k21*PERIPH;
$TABLE double CP = CENT/V1;
$CAPTURE CP
", compile = TRUE)

# 3-compartment IV
threecomp_iv <- mcode("threecomp_iv", "
$PARAM CL = 5, V1 = 50, Q2 = 8, V2 = 80, Q3 = 3, V3 = 150
$CMT CENT PERIPH1 PERIPH2
$ODE
double k10 = CL/V1;
double k12 = Q2/V1;
double k21 = Q2/V2;
double k13 = Q3/V1;
double k31 = Q3/V3;
dxdt_CENT = -k10*CENT - k12*CENT + k21*PERIPH1 - k13*CENT + k31*PERIPH2;
dxdt_PERIPH1 = k12*CENT - k21*PERIPH1;
dxdt_PERIPH2 = k13*CENT - k31*PERIPH2;
$TABLE double CP = CENT/V1;
$CAPTURE CP
", compile = TRUE)

# Transit absorption
transit_abs <- mcode("transit_abs", "
$PARAM N = 3, KTR = 0.5, KA = 1.5, CL = 5, V = 50
$CMT T1 T2 T3 ABS CENT
$ODE
dxdt_T1 = -KTR * T1;
dxdt_T2 = KTR * T1 - KTR * T2;
dxdt_T3 = KTR * T2 - KTR * T3;
dxdt_ABS = KTR * T3 - KA * ABS;
dxdt_CENT = KA * ABS - CL/V * CENT;
$TABLE double CP = CENT/V;
$CAPTURE CP
", compile = TRUE)

# Michaelis-Menten
mm_elim <- mcode("mm_elim", "
$PARAM VMAX = 50, KM = 10, V = 50
$CMT CENT
$ODE
double C = CENT/V;
dxdt_CENT = -VMAX * C / (KM + C);
$TABLE double CP = CENT/V;
$CAPTURE CP
", compile = TRUE)

# Direct Emax PD
pkpd_emax <- mcode("pkpd_emax", "
$PARAM CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10
$CMT CENT
$ODE dxdt_CENT = -CL/V * CENT;
$TABLE
double CP = CENT/V;
double EFF = E0 + EMAX * CP / (EC50 + CP);
$CAPTURE CP EFF
", compile = TRUE)

# Sigmoid Emax PD
pkpd_sigmoid <- mcode("pkpd_sigmoid", "
$PARAM CL = 5, V = 50, E0 = 0, EMAX = 100, EC50 = 10, GAMMA = 2
$CMT CENT
$ODE dxdt_CENT = -CL/V * CENT;
$TABLE
double CP = CENT/V;
double CG = pow(CP, GAMMA);
double EC50G = pow(EC50, GAMMA);
double EFF = E0 + EMAX * CG / (EC50G + CG);
$CAPTURE CP EFF
", compile = TRUE)

# Indirect response IRM1
pkpd_irm1 <- mcode("pkpd_irm1", "
$PARAM CL = 5, V = 50, KIN = 10, KOUT = 0.1, IMAX = 0.9, IC50 = 5
$CMT CENT RESPONSE
$MAIN RESPONSE_0 = KIN/KOUT;
$ODE
double CP = CENT/V;
double INH = IMAX * CP / (IC50 + CP);
dxdt_CENT = -CL/V * CENT;
dxdt_RESPONSE = KIN * (1 - INH) - KOUT * RESPONSE;
$CAPTURE CP
", compile = TRUE)

# Biophase equilibration
pkpd_biophase <- mcode("pkpd_biophase", "
$PARAM CL = 5, V = 50, KE0 = 0.5, E0 = 0, EMAX = 100, EC50 = 10
$CMT CENT EFFECT
$ODE
double CP = CENT/V;
dxdt_CENT = -CL/V * CENT;
dxdt_EFFECT = KE0 * (CP - EFFECT);
$TABLE
double CE = EFFECT;
double EFF = E0 + EMAX * CE / (EC50 + CE);
$CAPTURE CP CE EFF
", compile = TRUE)

# Population model with IIV
twocomp_pop <- mcode("twocomp_pop", "
$PARAM TVCL = 5, TVV1 = 50, TVQ = 10, TVV2 = 100
$PARAM ETA1 = 0, ETA2 = 0
$CMT CENT PERIPH
$MAIN
double CL = TVCL * exp(ETA1);
double V1 = TVV1 * exp(ETA2);
double Q = TVQ;
double V2 = TVV2;
$ODE
double k10 = CL/V1;
double k12 = Q/V1;
double k21 = Q/V2;
dxdt_CENT = -k10*CENT - k12*CENT + k21*PERIPH;
dxdt_PERIPH = k12*CENT - k21*PERIPH;
$TABLE double CP = CENT/V1;
$CAPTURE CP CL V1
", compile = TRUE)

cat("All models compiled successfully.\n\n")

# =============================================================================
# Benchmark Execution
# =============================================================================

all_results <- data.frame()

# -----------------------------------------------------------------------------
# 1. PK Models
# -----------------------------------------------------------------------------
cat(strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PK MODELS\n")
cat(strrep("=", 70), "\n")

# Common settings
dose_data <- data.frame(ID = 1, time = 0, amt = 100, cmt = 1, evid = 1)
times <- seq(0, 24, by = 0.5)

# One-compartment IV
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(onecomp_iv, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "compartmental", "OneCompIVBolus"
))

# One-compartment oral
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(onecomp_oral, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "compartmental", "OneCompOralFirstOrder"
))

# Two-compartment IV
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(twocomp_iv, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "compartmental", "TwoCompIVBolus"
))

# Two-compartment oral
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(twocomp_oral, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "compartmental", "TwoCompOral"
))

# Three-compartment IV
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(threecomp_iv, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "compartmental", "ThreeCompIVBolus"
))

# Transit absorption
dose_transit <- data.frame(ID = 1, time = 0, amt = 100, cmt = 1, evid = 1)
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(transit_abs, dose_transit, end = 24, delta = 0.5)),
  "simulate", "pk", "absorption", "TransitAbsorption"
))

# Michaelis-Menten
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(mm_elim, dose_data, end = 24, delta = 0.5)),
  "simulate", "pk", "nonlinear", "MichaelisMenten"
))

# Multiple dosing (7 days)
dose_multi <- data.frame(
  ID = 1,
  time = seq(0, 144, by = 24),
  amt = 100,
  cmt = 1,
  evid = 1
)
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(twocomp_iv, dose_multi, end = 168, delta = 2)),
  "simulate_multiple_dose", "pk", "dosing", "TwoCompIV_7days"
))

# -----------------------------------------------------------------------------
# 2. PD Models
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PD MODELS\n")
cat(strrep("=", 70), "\n")

# Direct Emax
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_emax, dose_data, end = 24, delta = 0.5)),
  "direct_effect", "pd", "direct", "DirectEmax"
))

# Sigmoid Emax
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_sigmoid, dose_data, end = 24, delta = 0.5)),
  "direct_effect", "pd", "direct", "SigmoidEmax"
))

# IRM1
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_irm1, dose_data, end = 72, delta = 1)),
  "indirect_response", "pd", "indirect", "IRM1"
))

# Biophase
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_biophase, dose_data, end = 72, delta = 1)),
  "effect_compartment", "pd", "biophase", "BiophaseEquilibration"
))

# -----------------------------------------------------------------------------
# 3. PKPD Coupled
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: PKPD COUPLED\n")
cat(strrep("=", 70), "\n")

# These are the same models but categorized as PKPD
all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_emax, dose_data, end = 72, delta = 1)),
  "pkpd_coupled", "pkpd", "direct", "OneComp_DirectEmax"
))

all_results <- rbind(all_results, run_benchmark(
  quote(mrgsim(pkpd_sigmoid, dose_data, end = 72, delta = 1)),
  "pkpd_coupled", "pkpd", "direct", "OneComp_SigmoidEmax"
))

# -----------------------------------------------------------------------------
# 4. NCA (computed manually from simulation output)
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: NCA\n")
cat(strrep("=", 70), "\n")

# Generate PK data for NCA
sim_result <- mrgsim(twocomp_iv, dose_data, end = 24, delta = 0.5)
pk_data <- as.data.frame(sim_result)
times_nca <- pk_data$time
conc_nca <- pk_data$CP

# Cmax/Tmax
all_results <- rbind(all_results, run_benchmark(
  quote({
    cmax <- max(conc_nca)
    tmax <- times_nca[which.max(conc_nca)]
    c(cmax, tmax)
  }),
  "cmax_tmax", "nca", "exposure", "Cmax_Tmax"
))

# AUC (trapezoidal)
all_results <- rbind(all_results, run_benchmark(
  quote({
    n <- length(times_nca)
    auc <- sum((conc_nca[-1] + conc_nca[-n]) / 2 * diff(times_nca))
    auc
  }),
  "auc_0_t", "nca", "exposure", "AUC_0_t"
))

# Lambda-z estimation (last 3+ points)
all_results <- rbind(all_results, run_benchmark(
  quote({
    # Use last 5 points for terminal slope
    n <- length(conc_nca)
    idx <- (n-4):n
    lm_fit <- lm(log(conc_nca[idx] + 1e-10) ~ times_nca[idx])
    lambda_z <- -coef(lm_fit)[2]
    lambda_z
  }),
  "lambda_z", "nca", "kinetics", "Lambda_z"
))

# Half-life
all_results <- rbind(all_results, run_benchmark(
  quote({
    n <- length(conc_nca)
    idx <- (n-4):n
    lm_fit <- lm(log(conc_nca[idx] + 1e-10) ~ times_nca[idx])
    lambda_z <- -coef(lm_fit)[2]
    t_half <- log(2) / lambda_z
    t_half
  }),
  "half_life", "nca", "kinetics", "Half_life"
))

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
  omega <- matrix(c(0.09, 0, 0, 0.04), 2, 2)  # 30% CV CL, 20% CV V
  eta <- MASS::mvrnorm(n_subj, c(0, 0), omega)

  pop_data <- data.frame(
    ID = rep(1:n_subj, each = 1),
    time = 0,
    amt = 100,
    cmt = 1,
    evid = 1,
    ETA1 = eta[, 1],
    ETA2 = eta[, 2]
  )

  all_results <- rbind(all_results, run_benchmark(
    quote(mrgsim(twocomp_pop, pop_data, end = 24, delta = 0.5)),
    "population_simulation", "population", "iiv",
    "TwoCompIV_IIV", n_subjects = n_subj
  ))
}

# -----------------------------------------------------------------------------
# 6. Sensitivity Analysis
# -----------------------------------------------------------------------------
cat("\n", strrep("=", 70), "\n")
cat("BENCHMARK CATEGORY: SENSITIVITY ANALYSIS\n")
cat(strrep("=", 70), "\n")

# Local sensitivity (perturbation)
all_results <- rbind(all_results, run_benchmark(
  quote({
    # Baseline
    base <- mrgsim(onecomp_iv, dose_data, end = 24, delta = 0.5)

    # Perturbed CL (+10%)
    pert_mod <- param(onecomp_iv, CL = 5.5)
    pert <- mrgsim(pert_mod, dose_data, end = 24, delta = 0.5)

    # Calculate sensitivity
    sens <- (as.data.frame(pert)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP
    sens
  }),
  "local_single", "sensitivity", "local", "LocalSingle"
))

# Multi-parameter sensitivity
all_results <- rbind(all_results, run_benchmark(
  quote({
    base <- mrgsim(onecomp_iv, dose_data, end = 24, delta = 0.5)

    # Perturb CL
    pert_cl <- mrgsim(param(onecomp_iv, CL = 5.5), dose_data, end = 24, delta = 0.5)
    sens_cl <- (as.data.frame(pert_cl)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP

    # Perturb V
    pert_v <- mrgsim(param(onecomp_iv, V = 55), dose_data, end = 24, delta = 0.5)
    sens_v <- (as.data.frame(pert_v)$CP - as.data.frame(base)$CP) / as.data.frame(base)$CP

    list(sens_cl, sens_v)
  }),
  "local_multi", "sensitivity", "local", "LocalMulti"
))

# =============================================================================
# Save Results
# =============================================================================

output_file <- file.path(CONFIG$output_dir, "mrgsolve_comprehensive_benchmarks.csv")
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
