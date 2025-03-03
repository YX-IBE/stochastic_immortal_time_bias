setwd("~/S1")

library(tidyverse)
library(data.table)
library(gfoRmula)
library(hash)
library(doParallel)
library(foreach)

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 2)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
data.table::setDTthreads(1)

# Set up parallel backend
n_cores <- 120
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
n_batch <- 10
batch_size <- 100

# Var1 and Var2
df_var <- as.data.frame(cbind(ID = 1 : (n_pat * 4),  Var1 = rep(NA, (n_pat * 4)), Var2 = rep(NA, (n_pat * 4))))
df_var$Var1 <- rep(rep(c(0, 1), 2), n_pat)
df_var$Var2 <- rep(rep(c(0, 1), c(2, 2)), n_pat)

# G-formula
# Sequential computation for each path_list in each batch 100 * 25
# Parallel bootstrap for each model fit

#######################
# Time scale: forward #
#######################

dirs <- c("./model_fit",
          "./model_fit/forward", 
          "./model_fit/forward/gformula", 
          "./model_fit/forward/gformula/input", 
          "./model_fit/forward/gformula/output")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

# Parameters for g-formula

source("./code/functions_Markov_continuous.R")

for (batch in 1 : n_batch) {
  path_lists_batch <- readRDS(paste0("./simulation/forward/path_lists_batch_", batch, ".rds"))
  cat("Prepare for batch", batch, "\n")
  
  foreach(simul = 1:batch_size, .packages = c("tidyverse", "data.table", "gfoRmula", "hash")) %dopar% {
    
    path_list <- path_lists_batch[[simul]]
    
    bootstrap_samples <- 1000
    seed <- 1809
    
    t_interval <- 1
    gform_dat <- create_gform_dat(path_list = path_list, t_star = t_star, t_interval = t_interval)
    gform_dat <- gform_dat %>% left_join(df_var, by = "ID")
    gform_dat <- data.table(gform_dat)
    
    id <- "ID"
    time_name <- "Time"
    time_points <- t_star / t_interval + 1
    outcome_name <- "Death"
    outcome_type <- "survival"
    basecovs <- c("Var1", "Var2")
    covnames <- c("Treat", "Relapse")
    covtypes <- c("absorbing", "absorbing")
    histories <- c(lagged, cumavg)
    histvars <- list(c("Treat", "Relapse"), c("Treat", "Relapse"))
    covparams <- list(covmodels = c(Treat ~ Time,
                                    Relapse ~ Treat + Var1 + lag1_Treat + cumavg_Treat + Time))
    restrictions <- list(c("Treat", "Relapse == 0", carry_forward))
    ymodel <- Death ~ Relapse + Treat + Var1 + Var2 + lag1_Relapse + lag1_Treat + cumavg_Relapse + cumavg_Treat + Time
    
    intvars <- rep(list("Treat"), 6)
    interventions <- NULL
    int_descript <- NULL
    
    custom_int <- function(newdf, pool, intvar, intvals, time_name, t){
      if (t == 0) {
        newdf[, rel_time := 0]
      } else {
        newdf[rel_time == 0 & Relapse == 1, rel_time := t]
      }
      newdf[, (intvar) := 0]
      newdf[(rel_time == 0 | rel_time > intvals[[1]]) & t >= intvals[[1]], (intvar) := 1]
    }
    
    
    for (i in c(1, 3, 6, 9, 12, 37)) {
      interventions <- c(interventions, list(list(c(custom_int, i))))
      if (i == 1) {
        int_descript <- c(int_descript, "treat immediately")
      } else if (i == 37) {
        int_descript <- c(int_descript, "never treat")
      } else {
        int_descript <- c(int_descript, paste("treat at", i, "month"))
      }
    }
    
    intcomp <- list(6, 1) # hazard ratio, 6 (never treat) as the reference intervention
    model_fits <- T
    ci_method <- "percentile"
    
    params <- hash()
    params[["obs_data"]] <- gform_dat
    params[["id"]] <- id
    params[["time_name"]] <- time_name
    params[["time_points"]] <- time_points
    params[["outcome_name"]] <- outcome_name
    params[["outcome_type"]] <- outcome_type
    params[["basecovs"]] <- basecovs
    params[["covnames"]] <- covnames
    params[["covtypes"]] <- covtypes
    params[["histories"]] <- histories
    params[["histvars"]] <- histvars
    params[["covparams"]] <- covparams
    params[["restrictions"]] <- restrictions
    params[["ymodel"]] <- ymodel
    params[["model_fits"]] <- model_fits
    params[["ci_method"]] <- ci_method
    params[["intcomp"]] <- intcomp
    params[["intvars"]] <- intvars
    params[["interventions"]] <- interventions
    params[["ref_int"]] <- length(int_descript)
    params[["int_descript"]] <- int_descript
    params[["seed"]] <- seed
    params[["nsamples"]] <- bootstrap_samples
    
    saveRDS(params, file = paste0("./model_fit/forward/gformula/input/gformula_batch_", batch, "_simul_", simul, ".rds"))
  }
}

#######################
# Time scale: renewal #
#######################

dirs <- c("./model_fit",
          "./model_fit/renewal", 
          "./model_fit/renewal/gformula", 
          "./model_fit/renewal/gformula/input", 
          "./model_fit/renewal/gformula/output")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

# Parameters for g-formula

source("./code/functions_Markov_continuous.R")

for (batch in 1 : n_batch) {
  path_lists_batch <- readRDS(paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"))
  cat("Prepare for batch", batch, "\n")
  
  foreach(simul = 1:batch_size, .packages = c("tidyverse", "data.table", "gfoRmula", "hash")) %dopar% {
    
    path_list <- path_lists_batch[[simul]]
    
    bootstrap_samples <- 1000
    seed <- 1809
    
    t_interval <- 1
    gform_dat <- create_gform_dat(path_list = path_list, t_star = t_star, t_interval = t_interval)
    gform_dat <- gform_dat %>% left_join(df_var, by = "ID")
    gform_dat <- data.table(gform_dat)
    
    id <- "ID"
    time_name <- "Time"
    time_points <- t_star / t_interval + 1
    outcome_name <- "Death"
    outcome_type <- "survival"
    basecovs <- c("Var1", "Var2")
    covnames <- c("Treat", "Relapse")
    covtypes <- c("absorbing", "absorbing")
    histories <- c(lagged, cumavg)
    histvars <- list(c("Treat", "Relapse"), c("Treat", "Relapse"))
    covparams <- list(covmodels = c(Treat ~ Time,
                                    Relapse ~ Treat + Var1 + lag1_Treat + cumavg_Treat + Time))
    restrictions <- list(c("Treat", "Relapse == 0", carry_forward))
    ymodel <- Death ~ Relapse + Treat + Var1 + Var2 + lag1_Relapse + lag1_Treat + cumavg_Relapse + cumavg_Treat + Time
    
    intvars <- rep(list("Treat"), 6)
    interventions <- NULL
    int_descript <- NULL
    
    custom_int <- function(newdf, pool, intvar, intvals, time_name, t){
      if (t == 0) {
        newdf[, rel_time := 0]
      } else {
        newdf[rel_time == 0 & Relapse == 1, rel_time := t]
      }
      newdf[, (intvar) := 0]
      newdf[(rel_time == 0 | rel_time > intvals[[1]]) & t >= intvals[[1]], (intvar) := 1]
    }
    
    
    for (i in c(1, 3, 6, 9, 12, 37)) {
      interventions <- c(interventions, list(list(c(custom_int, i))))
      if (i == 1) {
        int_descript <- c(int_descript, "treat immediately")
      } else if (i == 37) {
        int_descript <- c(int_descript, "never treat")
      } else {
        int_descript <- c(int_descript, paste("treat at", i, "month"))
      }
    }
    
    intcomp <- list(6, 1) # hazard ratio, 6 (never treat) as the reference intervention
    model_fits <- T
    ci_method <- "percentile"
    
    params <- hash()
    params[["obs_data"]] <- gform_dat
    params[["id"]] <- id
    params[["time_name"]] <- time_name
    params[["time_points"]] <- time_points
    params[["outcome_name"]] <- outcome_name
    params[["outcome_type"]] <- outcome_type
    params[["basecovs"]] <- basecovs
    params[["covnames"]] <- covnames
    params[["covtypes"]] <- covtypes
    params[["histories"]] <- histories
    params[["histvars"]] <- histvars
    params[["covparams"]] <- covparams
    params[["restrictions"]] <- restrictions
    params[["ymodel"]] <- ymodel
    params[["model_fits"]] <- model_fits
    params[["ci_method"]] <- ci_method
    params[["intcomp"]] <- intcomp
    params[["intvars"]] <- intvars
    params[["interventions"]] <- interventions
    params[["ref_int"]] <- length(int_descript)
    params[["int_descript"]] <- int_descript
    params[["seed"]] <- seed
    params[["nsamples"]] <- bootstrap_samples
    
    saveRDS(params, file = paste0("./model_fit/renewal/gformula/input/gformula_batch_", batch, "_simul_", simul, ".rds"))
  }
}

#######################
# Time scale: extended #
#######################

dirs <- c("./model_fit",
          "./model_fit/extended", 
          "./model_fit/extended/gformula", 
          "./model_fit/extended/gformula/input", 
          "./model_fit/extended/gformula/output")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

# Parameters for g-formula

source("./code/functions_Markov_continuous.R")

for (batch in 1 : n_batch) {
  path_lists_batch <- readRDS(paste0("./simulation/extended/path_lists_batch_", batch, ".rds"))
  cat("Prepare for batch", batch, "\n")
  
  foreach(simul = 1:batch_size, .packages = c("tidyverse", "data.table", "gfoRmula", "hash")) %dopar% {
    
    path_list <- path_lists_batch[[simul]]
    
    bootstrap_samples <- 1000
    seed <- 1809
    
    t_interval <- 1
    gform_dat <- create_gform_dat(path_list = path_list, t_star = t_star, t_interval = t_interval)
    gform_dat <- gform_dat %>% left_join(df_var, by = "ID")
    gform_dat <- data.table(gform_dat)
    
    id <- "ID"
    time_name <- "Time"
    time_points <- t_star / t_interval + 1
    outcome_name <- "Death"
    outcome_type <- "survival"
    basecovs <- c("Var1", "Var2")
    covnames <- c("Treat", "Relapse")
    covtypes <- c("absorbing", "absorbing")
    histories <- c(lagged, cumavg)
    histvars <- list(c("Treat", "Relapse"), c("Treat", "Relapse"))
    covparams <- list(covmodels = c(Treat ~ Time,
                                    Relapse ~ Treat + Var1 + lag1_Treat + cumavg_Treat + Time))
    restrictions <- list(c("Treat", "Relapse == 0", carry_forward))
    ymodel <- Death ~ Relapse + Treat + Var1 + Var2 + lag1_Relapse + lag1_Treat + cumavg_Relapse + cumavg_Treat + Time
    
    intvars <- rep(list("Treat"), 6)
    interventions <- NULL
    int_descript <- NULL
    
    custom_int <- function(newdf, pool, intvar, intvals, time_name, t){
      if (t == 0) {
        newdf[, rel_time := 0]
      } else {
        newdf[rel_time == 0 & Relapse == 1, rel_time := t]
      }
      newdf[, (intvar) := 0]
      newdf[(rel_time == 0 | rel_time > intvals[[1]]) & t >= intvals[[1]], (intvar) := 1]
    }
    
    
    for (i in c(1, 3, 6, 9, 12, 37)) {
      interventions <- c(interventions, list(list(c(custom_int, i))))
      if (i == 1) {
        int_descript <- c(int_descript, "treat immediately")
      } else if (i == 37) {
        int_descript <- c(int_descript, "never treat")
      } else {
        int_descript <- c(int_descript, paste("treat at", i, "month"))
      }
    }
    
    intcomp <- list(6, 1) # hazard ratio, 6 (never treat) as the reference intervention
    model_fits <- T
    ci_method <- "percentile"
    
    params <- hash()
    params[["obs_data"]] <- gform_dat
    params[["id"]] <- id
    params[["time_name"]] <- time_name
    params[["time_points"]] <- time_points
    params[["outcome_name"]] <- outcome_name
    params[["outcome_type"]] <- outcome_type
    params[["basecovs"]] <- basecovs
    params[["covnames"]] <- covnames
    params[["covtypes"]] <- covtypes
    params[["histories"]] <- histories
    params[["histvars"]] <- histvars
    params[["covparams"]] <- covparams
    params[["restrictions"]] <- restrictions
    params[["ymodel"]] <- ymodel
    params[["model_fits"]] <- model_fits
    params[["ci_method"]] <- ci_method
    params[["intcomp"]] <- intcomp
    params[["intvars"]] <- intvars
    params[["interventions"]] <- interventions
    params[["ref_int"]] <- length(int_descript)
    params[["int_descript"]] <- int_descript
    params[["seed"]] <- seed
    params[["nsamples"]] <- bootstrap_samples
    
    saveRDS(params, file = paste0("./model_fit/extended/gformula/input/gformula_batch_", batch, "_simul_", simul, ".rds"))
  }
}

stopCluster(cl)
registerDoSEQ()


