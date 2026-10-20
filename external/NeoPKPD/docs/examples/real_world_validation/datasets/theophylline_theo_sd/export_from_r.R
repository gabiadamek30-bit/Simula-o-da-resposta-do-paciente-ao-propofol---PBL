# Export theo_sd from nlmixr2data to CSV
# This script is intentionally boring and explicit.

if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  install.packages("nlmixr2data", repos = "https://cloud.r-project.org")
}

library(nlmixr2data)

d <- nlmixr2data::theo_sd

# Enforce expected columns and order
cols <- c("ID","TIME","DV","AMT","EVID","CMT","WT")
missing <- setdiff(cols, names(d))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

d <- d[, cols]

out <- "theo_sd.csv"
write.csv(d, out, row.names = FALSE, quote = TRUE)
cat("Wrote", out, "\n")
