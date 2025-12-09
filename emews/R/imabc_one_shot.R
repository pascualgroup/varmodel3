library(data.table)
.datatable.aware = TRUE
assignInNamespace("cedta.pkgEvalsUserCode", c(data.table:::cedta.pkgEvalsUserCode,"imabc"), "data.table")

library(jsonlite)
library(lhs)
library(imabc)
library(doParallel)
library(truncnorm)
library(readr)
library(dplyr)
library(yaml)

options(stringsAsFactors=F)
options(error=function()traceback(2))
options(imabc.target_eval_distance = "zscore")
# options(datatable.verbose = TRUE)


string_to_list_of_vectors <- function(x){
  lapply(unlist(strsplit(x,";")),function(y) as.numeric(unlist(strsplit(y,","))))
}

instance <- 1

b_fun <- function(params, all_parm_names, target_fun, other_inputs) {
  nrows <- nrow(params)
  # make a copy so that we can add / remove columns
  # without causing issues in imabc later
  parms_to_run <- copy(params)
  current.iter <- parms_to_run[, iter][1]
  print(current.iter)
  parms_to_run[, c("instance") := list(seq(instance, length.out=nrows))]
  saveRDS(parms_to_run, file=paste0(output_dir, "/parms_to_run_", current.iter, ".RDS"))
  parms_to_run[, c("iter", "draw", "step", "scaled_dist", "sample_wt") := list(NULL, NULL, NULL, NULL, NULL)]
  instance <<- instance + nrows
  
  # print(parms_to_run)
  # json_params <- toJSON(parms_to_run, auto_unbox = T, digits = NA, rownames = F)
  # json_lower_bounds <- toJSON(as.list(other_inputs$targets$current_lower_bounds), digits = NA, auto_unbox = T)
  # json_upper_bounds <- toJSON(as.list(other_inputs$targets$current_upper_bounds), digits = NA, auto_unbox = T)
  # cat(paste0("Out Start ", format(Sys.time(), "%FT%T%z"), "\n"))
  # print(paste0("json: ", json_params))
  # Push string representation of parms.to.run to queue
  # OUT_put(json_params)

  # ';' separated list of json maps
  # res <- IN_get()
  # res_list <- lapply(unlist(strsplit(res, ";")), function(x) fromJSON(x))
  # Ynew <- do.call(rbind, res_list)
  # saveRDS(Ynew, file=paste0(turbine_output, "/Ynew_", current.iter, ".RDS"))
  # return(as.data.frame(Ynew))
  if (current.iter == end_iter) {
    print("DONE")
    return()
  } else {
    results_file <- cfg$results_files[current.iter]
    Ynew <- read.csv(results_file)
    Ynew <- Ynew %>% arrange(instance)
    Ynew <- Ynew[, !(names(Ynew) %in% "instance")]
    print(head(Ynew))
    return(Ynew)
  }
}

cat(paste0("IMABC Start ", format(Sys.time(), "%FT%T%z"), "\n"))

args <- commandArgs(trailingOnly = TRUE)
cfg_file = args[1]
cfg = read_yaml(cfg_file)
source(cfg$algo_param_file)

# priors are placed here by *run_imabc.sh
priors.path = cfg$priors
priors.df <- data.frame(read_csv(priors.path))
priors <- as.priors(priors.df)

targets.path <- cfg$targets
raw_df <- data.frame(read_csv(targets.path))
if (cfg$n_targets == 4) {
  # "prevalence","meanMOIvar","meanPTS","inverseSimpsonIndex"
  targets_df = subset(raw_df, target_names %in% c("prevalence","meanMOIvar","meanPTS","inverseSimpsonIndex"))
} else if (cfg$n_targets == 6) {
  targets_df = subset(raw_df, target_names %in% c("prevalence","meanMOIvar","meanPTS","meanPTSGroupBC","inverseSimpsonIndex","inverseSimpsonIndexBC"))
} else {
  stop("Bad number of targets")
}
targets <- as.targets(targets_df)
print(targets)

print(cfg$results_files)


target_suffix <- tools::file_path_sans_ext(basename(targets.path))
# output_dir <- paste0(cfg$output_directory, "_", target_suffix)
# dir.create(output_dir, FALSE)

output_dir <- cfg$output_directory
setwd(output_dir)
end_iter <- cfg$end_iter

algo.params$imabc.args$max_iter = cfg$end_iter
# use modifyList to override items in algo.params$imabc.args
imabc.args <- modifyList(algo.params$imabc.args, list(output_directory = output_dir,
                         targets=targets, priors = priors, backend_fun = b_fun))

# print(imabc.args$priors)
# print(imabc.args$targets)

a <- do.call(imabc, imabc.args)

