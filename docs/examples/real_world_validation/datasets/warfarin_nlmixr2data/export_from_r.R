if (!requireNamespace("nlmixr2data", quietly = TRUE)) {
  install.packages("nlmixr2data", repos = "https://cloud.r-project.org")
}
library(nlmixr2data)

d <- nlmixr2data::warfarin

cols <- c("id","time","amt","dv","dvid","evid","wt","age","sex")
missing <- setdiff(cols, names(d))
if (length(missing) > 0) stop(paste("Missing columns:", paste(missing, collapse=", ")))

d <- d[, cols]

out <- "warfarin.csv"
write.csv(d, out, row.names = FALSE, quote = TRUE)
cat("Wrote", out, "\n")
