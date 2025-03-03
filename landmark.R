# Landmark supermodel
setwd("~/S1")

library(foreach)
library(doParallel)
library(tidyverse)
library(survival)
library(mstate)

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

n_cores <- 100
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
n_batch <- 10
batch_size <- 100
scales <- c("forward", "renewal", "extended")
int_list <- c(1, 3, 6, 9, 12)
sim_pat <- 1000

# Var1 and Var2
df_var <- as.data.frame(cbind(ID = 1 : (n_pat * 4),  Var1 = rep(NA, (n_pat * 4)), Var2 = rep(NA, (n_pat * 4))))
df_var$Var1 <- rep(rep(c(0, 1), 2), n_pat)
df_var$Var2 <- rep(rep(c(0, 1), c(2, 2)), n_pat)

# Parameters for super landmark data
seq_horizon <- c(6, 12)
seq_interval <- c(1, 3)

for (scale in scales) {
  ms_scale <- ifelse(scale == "extended", "renewal", scale)
  
  dirs <- c("./model_fit",
            paste0("./model_fit/", scale), 
            paste0("./model_fit/", scale, "/landmark"),
            paste0("./model_fit/", scale, "/landmark/int_1"),
            paste0("./model_fit/", scale, "/landmark/int_3"))
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir)
      print(paste("Directory", dir, "created."))
    } else {
      print(paste("Directory", dir, "already exists."))
    }
  }
  
  for (t_interval in seq_interval) {
    for (t_horizon in seq_horizon) {
      
      # Load each batch sequentially
      for (batch in 1 : n_batch) {
        
        file_name <- paste0("./model_fit/", scale, "/landmark/int_", t_interval, "/", ms_scale, "_hor_", t_horizon, "_batch_", batch, ".rds")
        
        # Skip processing if the file already exists
        if (file.exists(file_name)) {
          cat(sprintf("File already exists: %s. Skipping...\n", file_name))
          next
        }
        
        path_lists_batch <- readRDS(paste0("./simulation/", scale, "/path_lists_batch_", batch, ".rds"))
        
        # fit the model in parallel for each path_list in the current batch
        landmark_results_batch <- foreach(simul = 1 : batch_size, .combine = function(x, y) c(x, list(y)), .init = list(), .packages = c("tidyverse", "survival")) %dopar% {
          source("functions.R")
          
          path_list <- path_lists_batch[[simul]]
          
          sim_dat_int_0 <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_36_", ms_scale, ".rds"))
          sim_dat_int_0$Treat_s <- 0
          dat_landmark_int_0 <- create_landmark_int(sim_dat_int = sim_dat_int_0, int_time = t_star, t_star = t_star)
          col_ind <- which(colnames(dat_landmark_int_0) %in% c("Treat", "Relapse", "Var1", "Var2"))
          dat_landmark_int_0[, col_ind] <- lapply(dat_landmark_int_0[, col_ind], function(x) factor(x, levels = c(0, 1)))
          
          dat_landmark_int_n_list <- NULL
          for (int_time in int_list) {
            sim_dat_int_n <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_", int_time, "_", ms_scale, ".rds"))
            dat_landmark_int_n <- create_landmark_int(sim_dat_int = sim_dat_int_n, int_time = int_time, t_star = t_star)
            col_ind <- which(colnames(dat_landmark_int_n) %in% c("Treat", "Relapse", "Var1", "Var2"))
            dat_landmark_int_n[, col_ind] <- lapply(dat_landmark_int_n[, col_ind], function(x) factor(x, levels = c(0, 1)))
            dat_landmark_int_n_list <- c(dat_landmark_int_n_list, list(dat_landmark_int_n))
          }
          
          dat_landmark <- create_landmark_dat(path_list = path_list, t_star = t_star, t_interval = t_interval, t_horizon = t_horizon) 
          dat_landmark <- dat_landmark %>% left_join(df_var, by = "ID")
          
          col_ind <- which(colnames(dat_landmark) %in% c("Treat", "Relapse", "Var1", "Var2"))
          dat_landmark[, col_ind] <- lapply(dat_landmark[, col_ind], function(x) factor(x, levels = c(0, 1)))
          
          # Avoid zero-length interval between LM and Death_t due to high precision data
          threshold <- 1e-06
          for (i in seq_len(nrow(dat_landmark))) {
            if(abs(dat_landmark$Death_t[i] - dat_landmark$LM[i]) < threshold) {
              ind_row <- which(dat_landmark$ID == dat_landmark$ID[i] & dat_landmark$Death_t == dat_landmark$Death_t[i])
              dat_landmark$Death_t[ind_row] <- dat_landmark$LM[i] + threshold
            }
          }
          
          landmark_fit <- coxph(Surv(LM, Death_t, Death_s) ~ Treat + Treat:Landmark + Treat:I(Landmark^2) +
                                  Relapse + Relapse:Landmark + Relapse:I(Landmark^2) +
                                  Landmark + I(Landmark^2) + I(Landmark^3) + cluster(ID), data = dat_landmark, method = "breslow")
          
          # 1. model coefficients / not for causal interpretation
          landmark_coef <- summary(landmark_fit)$coefficients
          
          # 2. ratio of cumulative hazard
          # 2.1 baseline cumulative hazard function
          baseline_dat <- data.frame(Treat = factor(0, levels = c(0, 1)), 
                                     Relapse = factor(0, levels = c(0, 1)), 
                                     Var1 = factor(0, levels = c(0, 1)), 
                                     Var2 = factor(0, levels = c(0, 1)), 
                                     Landmark = 0)
          
          base_haz <- basehaz(landmark_fit, newdata = baseline_dat)
          base_haz <- rbind(c(0, 0), base_haz)
          haz_fun <- splinefun(base_haz$time, base_haz$hazard, method = "monoH.FC")
          
          # 2.2 calculate ratio of cumulative hazards, in reference to never treat
          # 2.2.1 never treat (reference group)
          dat_landmark_int_0$exp_lp <- predict(landmark_fit, type = "risk", newdata = dat_landmark_int_0)
          
          # calculate step-wise cumulative hazard
          dat_landmark_int_0$discret_haz <- haz_fun(dat_landmark_int_0$Landmark + 1) - haz_fun(dat_landmark_int_0$Landmark)
          dat_landmark_int_0$cumhaz_step <- dat_landmark_int_0$discret_haz * dat_landmark_int_0$exp_lp
          
          cumhaz_0 <- unlist(tapply(dat_landmark_int_0$cumhaz_step, dat_landmark_int_0$ID, FUN = sum))
          
          est_abs_risk <- mean(1 - exp(-cumhaz_0))
          
          # 2.2.2 treat under specific intervention strategies
          est_cumhaz_contrast <- numeric()
          est_risk_contrast <- numeric()
          
          for (int_time in int_list) {
            dat_landmark_int_1 <- dat_landmark_int_n_list[[which(int_list == int_time)]]
            dat_landmark_int_1$exp_lp <- predict(landmark_fit, type = "risk", newdata = dat_landmark_int_1)
            
            # calculate step-wise cumulative hazard
            dat_landmark_int_1$discret_haz <- haz_fun(dat_landmark_int_1$Landmark + 1) - haz_fun(dat_landmark_int_1$Landmark)
            dat_landmark_int_1$cumhaz_step <- dat_landmark_int_1$discret_haz * dat_landmark_int_1$exp_lp
            
            cumhaz_1 <- unlist(tapply(dat_landmark_int_1$cumhaz_step, dat_landmark_int_1$ID, FUN = sum))
            
            est_abs_risk <- c(est_abs_risk, mean(1 - exp(-cumhaz_1)))
            est_cumhaz_contrast[which(int_list == int_time)] <- mean(cumhaz_1)/mean(cumhaz_0)
            est_risk_contrast[which(int_list == int_time)] <- mean(1 - exp(-cumhaz_1))/mean(1 - exp(-cumhaz_0))
          }
          
          landmark_results <- list(landmark_coef = landmark_coef, 
                                   est_abs_risk = est_abs_risk,
                                   est_cumhaz_contrast = est_cumhaz_contrast,
                                   est_risk_contrast = est_risk_contrast)
          
          return(landmark_results)
        }
        
        file_name <- paste0("./model_fit/", scale, "/landmark/int_", t_interval, "/", ms_scale, "_hor_", t_horizon, "_batch_", batch, ".rds")
        saveRDS(landmark_results_batch, file = file_name)
      }
    }
  }
}

stopCluster(cl)
