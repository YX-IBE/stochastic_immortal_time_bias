setwd("~/S1")

library(gfoRmula)
library(hash)
library(parallel)
library(doParallel)
library(foreach)

#######################
# Time scale: forward #
#######################

# Check unprocessed tasks
n_batch <- 10
batch_size <- 100

input_directory <- "./model_fit/forward/gformula/input"
output_directory <- "./model_fit/forward/gformula/output"

input_files <- list.files(path = input_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)
output_files <- list.files(path = output_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)

unprocessed_data <- list()
for (batch in 1 : n_batch) {
  for (simul in 1 : batch_size) {
    input_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    output_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    if (!any(input_files == input_file)) {
      stop(paste("Input", input_file, "doesn't exist"))
    }
    if (any(output_files == output_file)) {
      next
    }
    unprocessed_data[[length(unprocessed_data) + 1]] <- list(input_file = input_file, output_file = output_file)
  }
}
print(paste("Found", length(unprocessed_data), "unprocessed inputs..."))

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 2)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
data.table::setDTthreads(1)

n_cores <- 112
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

foreach(data = unprocessed_data, .errorhandling = "stop") %dopar% {
  library(gfoRmula)
  
  input_path <- file.path(input_directory, data$input_file)
  output_path <- file.path(output_directory, data$output_file)
  
  params <- readRDS(input_path)
  
  print(system.time({
    fit_gform <- gformula(obs_data = params[["obs_data"]], id = params[["id"]],
                          time_name = params[["time_name"]], time_points = params[["time_points"]],
                          outcome_name = params[["outcome_name"]], outcome_type = params[["outcome_type"]],
                          basecovs = params[["basecovs"]],
                          covnames = params[["covnames"]], covtypes = params[["covtypes"]], 
                          histories = params[["histories"]], histvars = params[["histvars"]], 
                          covparams = params[["covparams"]], restrictions = params[["restrictions"]],
                          ymodel = params[["ymodel"]], model_fits = params[["model_fits"]], ci_method = params[["ci_method"]],
                          intcomp = params[["intcomp"]], intvars = params[["intvars"]], interventions = params[["interventions"]], 
                          ref_int = params[["ref_int"]], int_descript = params[["int_descript"]], 
                          seed = params[["seed"]], nsamples = params[["nsamples"]], 
                          parallel = FALSE)
    saveRDS(fit_gform, file = output_path)
    }))
}

#######################
# Time scale: renewal #
#######################

# Check unprocessed tasks
n_batch <- 10
batch_size <- 100

input_directory <- "./model_fit/renewal/gformula/input"
output_directory <- "./model_fit/renewal/gformula/output"

input_files <- list.files(path = input_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)
output_files <- list.files(path = output_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)

unprocessed_data <- list()
for (batch in 1 : n_batch) {
  for (simul in 1 : batch_size) {
    input_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    output_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    if (!any(input_files == input_file)) {
      stop(paste("Input", input_file, "doesn't exist"))
    }
    if (any(output_files == output_file)) {
      next
    }
    unprocessed_data[[length(unprocessed_data) + 1]] <- list(input_file = input_file, output_file = output_file)
  }
}
print(paste("Found", length(unprocessed_data), "unprocessed inputs..."))

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 2)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
data.table::setDTthreads(1)

n_cores <- 26
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

foreach(data = unprocessed_data, .errorhandling = "stop") %dopar% {
  library(gfoRmula)
  
  input_path <- file.path(input_directory, data$input_file)
  output_path <- file.path(output_directory, data$output_file)
  
  params <- readRDS(input_path)
  
  print(system.time({
    fit_gform <- gformula(obs_data = params[["obs_data"]], id = params[["id"]],
                          time_name = params[["time_name"]], time_points = params[["time_points"]],
                          outcome_name = params[["outcome_name"]], outcome_type = params[["outcome_type"]],
                          basecovs = params[["basecovs"]],
                          covnames = params[["covnames"]], covtypes = params[["covtypes"]], 
                          histories = params[["histories"]], histvars = params[["histvars"]], 
                          covparams = params[["covparams"]], restrictions = params[["restrictions"]],
                          ymodel = params[["ymodel"]], model_fits = params[["model_fits"]], ci_method = params[["ci_method"]],
                          intcomp = params[["intcomp"]], intvars = params[["intvars"]], interventions = params[["interventions"]], 
                          ref_int = params[["ref_int"]], int_descript = params[["int_descript"]], 
                          seed = params[["seed"]], nsamples = params[["nsamples"]], 
                          parallel = FALSE)
    saveRDS(fit_gform, file = output_path)
  }))
}

#######################
# Time scale: extended #
#######################

# Check unprocessed tasks
n_batch <- 10
batch_size <- 100

input_directory <- "./model_fit/extended/gformula/input"
output_directory <- "./model_fit/extended/gformula/output"

input_files <- list.files(path = input_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)
output_files <- list.files(path = output_directory, pattern = NULL, all.files = FALSE, full.names = FALSE)

unprocessed_data <- list()
for (batch in 1 : n_batch) {
  for (simul in 1 : batch_size) {
    input_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    output_file <- paste0("gformula_batch_", batch, "_simul_", simul, ".rds")
    if (!any(input_files == input_file)) {
      stop(paste("Input", input_file, "doesn't exist"))
    }
    if (any(output_files == output_file)) {
      next
    }
    unprocessed_data[[length(unprocessed_data) + 1]] <- list(input_file = input_file, output_file = output_file)
  }
}
print(paste("Found", length(unprocessed_data), "unprocessed inputs..."))

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 2)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
data.table::setDTthreads(1)

n_cores <- 26
cluster <- makeCluster(n_cores)
registerDoParallel(cluster)

foreach(data = unprocessed_data, .errorhandling = "stop") %dopar% {
  library(gfoRmula)
  
  input_path <- file.path(input_directory, data$input_file)
  output_path <- file.path(output_directory, data$output_file)
  
  params <- readRDS(input_path)
  
  print(system.time({
    fit_gform <- gformula(obs_data = params[["obs_data"]], id = params[["id"]],
                          time_name = params[["time_name"]], time_points = params[["time_points"]],
                          outcome_name = params[["outcome_name"]], outcome_type = params[["outcome_type"]],
                          basecovs = params[["basecovs"]],
                          covnames = params[["covnames"]], covtypes = params[["covtypes"]], 
                          histories = params[["histories"]], histvars = params[["histvars"]], 
                          covparams = params[["covparams"]], restrictions = params[["restrictions"]],
                          ymodel = params[["ymodel"]], model_fits = params[["model_fits"]], ci_method = params[["ci_method"]],
                          intcomp = params[["intcomp"]], intvars = params[["intvars"]], interventions = params[["interventions"]], 
                          ref_int = params[["ref_int"]], int_descript = params[["int_descript"]], 
                          seed = params[["seed"]], nsamples = params[["nsamples"]], 
                          parallel = FALSE)
    saveRDS(fit_gform, file = output_path)
  }))
}

stopCluster(cluster)
stopImplicitCluster()
registerDoSEQ()


