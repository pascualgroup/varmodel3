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

options(stringsAsFactors=F)
options(error=function()traceback(2))
options(imabc.target_eval_distance = "zscore")
# options(datatable.verbose = TRUE)


emews_root <- Sys.getenv("EMEWS_PROJECT_ROOT")
if (emews_root == "") {
  r_root <- getwd()
} else {
  r_root <- paste0(emews_root, "/R")
}

turbine_output <- Sys.getenv("TURBINE_OUTPUT")
wd <- getwd()
setwd(turbine_output)

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
  parms_to_run[, c("instance") := list(seq(instance, length.out=nrows))]
  saveRDS(parms_to_run, file=paste0(turbine_output, "/parms_to_run_", current.iter, ".RDS"))
  parms_to_run[, c("iter", "draw", "step", "scaled_dist", "sample_wt") := list(NULL, NULL, NULL, NULL, NULL)]
  instance <<- instance + nrows
  
  # print(parms_to_run)
  json_params <- toJSON(parms_to_run, auto_unbox = T, digits = NA, rownames = F)
  # json_lower_bounds <- toJSON(as.list(other_inputs$targets$current_lower_bounds), digits = NA, auto_unbox = T)
  # json_upper_bounds <- toJSON(as.list(other_inputs$targets$current_upper_bounds), digits = NA, auto_unbox = T)
  # cat(paste0("Out Start ", format(Sys.time(), "%FT%T%z"), "\n"))
  # print(paste0("json: ", json_params))
  # Push string representation of parms.to.run to queue
  OUT_put(json_params)

  # ';' separated list of json maps
  res <- IN_get()
  res_list <- lapply(unlist(strsplit(res, ";")), function(x) fromJSON(x))
  Ynew <- do.call(rbind, res_list)
  saveRDS(Ynew, file=paste0(turbine_output, "/Ynew_", current.iter, ".RDS"))
  return(as.data.frame(Ynew))
}

cat(paste0("IMABC Start ", format(Sys.time(), "%FT%T%z"), "\n"))

# ask for parameters from queue
OUT_put("Params")
# accepts arguments to main_function, e.g., "pp = 2, it = 5"
res <- IN_get()

# TODO: source an R file, passed in through res "algo_params.R" instead of in the *.swift file
algo.params.file <- res
source(algo.params.file)

# priors are placed here by *run_imabc.sh
priors.path = paste0(turbine_output,"/priors.csv")
priors.df <- data.frame(read_csv(priors.path))
priors <- as.priors(priors.df)

targets.path <- paste0(turbine_output,"/targets.csv")
targets_df <- data.frame(read_csv(targets.path))
targets <- as.targets(targets_df)

# use modifyList to override items in algo.params$imabc.args
imabc.args <- modifyList(algo.params$imabc.args, list(output_directory = turbine_output,
                         targets=targets, priors = priors, backend_fun = b_fun))

# print(imabc.args$priors)
# print(imabc.args$targets)

a <- do.call(imabc,imabc.args)

OUT_put("DONE")

# This will be saved to experiment directory
saveRDS(a, file=paste0(turbine_output, "/imabc_res.Rds"))

OUT_put("Look at imabc_res.Rds for final results.")
cat(paste0("IMABC End ", format(Sys.time(), "%FT%T%z"), "\n"))
