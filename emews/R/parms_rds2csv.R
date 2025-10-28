library(data.table)

args <- commandArgs(trailingOnly = TRUE)
parms_file <- args[1]
parent <- dirname(parms_file)

dt <- readRDS(parms_file)
dt <- as.data.table(dt)
dt[, c("iter", "draw", "step", "scaled_dist", "sample_wt") := list(NULL, NULL, NULL, NULL, NULL)]
fwrite(dt, paste0(parent, "/parms_to_run.csv"))