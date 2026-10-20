if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  install.packages("nlmixr2data", repos = "https://cloud.r-project.org")
}
library(nlmixr2data)

d <- nlmixr2data::theo_md

cols <- c("ID","TIME","DV","AMT","EVID","CMT","WT")
missing <- setdiff(cols, names(d))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

d <- d[, cols]

out <- "theo_md.csv"
write.csv(d, out, row.names = FALSE, quote = TRUE)
cat("Wrote", out, "\n")
