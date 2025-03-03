# Time-dependent Cox

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

for (scale in scales) {
  ms_scale <- ifelse(scale == "extended", "renewal", scale)
  
  dirs <- c("./model_fit",
            paste0("./model_fit/", scale), 
            paste0("./model_fit/", scale, "/cox_td"))
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      print(paste("Directory", dir, "created."))
    } else {
      print(paste("Directory", dir, "already exists."))
    }
  }
  
  # model fit
  cox_fit_list_scale <- list(list("Surv(Start, Stop, Death) ~ Treat + Relapse + Var1 + Var2"), 
                             list("Surv(Start, Stop, Death) ~ Treat + Relapse + Var1 + Var2 + Start",
                                  "Surv(Start, Stop, Death) ~ Treat + Relapse + Var1 + Var2 + Start + I(Start^2)",
                                  "Surv(Start, Stop, Death) ~ Treat + Relapse + Var1 + Var2 + pspline(Start, df = 6)"))
  
  if (scale == "forward") {
    cox_fit_list <- cox_fit_list_scale[[1]]
  } else {
    cox_fit_list <- cox_fit_list_scale[[2]]
  }
  
  for (model in seq_len(length(cox_fit_list))) {
    
    # Load each batch sequentially  
    for (batch in 1 : n_batch) {
      
      if (scale == "forward") {
        file_name <- paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_fit_batch_", batch, ".rds")
      } else {
        file_name <- switch(model,
                            "1" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_linear_batch_", batch, ".rds"),
                            "2" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_quadratic_batch_", batch, ".rds"),
                            "3" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_pspline_batch_", batch, ".rds"))
      }
      
      # Skip processing if the file already exists
      if (file.exists(file_name)) {
        cat(sprintf("File already exists: %s. Skipping...\n", file_name))
        next
      }
      
      path_lists_batch <- readRDS(paste0("./simulation/", scale, "/path_lists_batch_", batch, ".rds"))
      
      # Fit the model in parallel for each path_list in the current batch    
      cox_results_batch <- foreach(simul = 1 : batch_size, .combine = function(x, y) c(x, list(y)), .init = list(), .packages = c("tidyverse", "survival")) %dopar% {
        source("functions.R")
        
        path_list <- path_lists_batch[[simul]]
        
        sim_dat_int_0 <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_36_", ms_scale, ".rds"))
        sim_dat_int_0$Treat_s <- 0
        dat_cox_int_0 <- create_cox_int(sim_dat_int = sim_dat_int_0, int_time = t_star, t_star = t_star)
        col_ind <- which(colnames(dat_cox_int_0) %in% c("Treat", "Relapse", "Var1", "Var2"))
        dat_cox_int_0[, col_ind] <- lapply(dat_cox_int_0[, col_ind], function(x) factor(x, levels = c(0, 1)))
        
        dat_cox_int_n_list <- NULL
        for (int_time in int_list) {
          sim_dat_int_n <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_", int_time, "_", ms_scale, ".rds"))
          dat_cox_int_n <- create_cox_int(sim_dat_int = sim_dat_int_n, int_time = int_time, t_star = t_star)
          col_ind <- which(colnames(dat_cox_int_n) %in% c("Treat", "Relapse", "Var1", "Var2"))
          dat_cox_int_n[, col_ind] <- lapply(dat_cox_int_n[, col_ind], function(x) factor(x, levels = c(0, 1)))
          dat_cox_int_n_list <- c(dat_cox_int_n_list, list(dat_cox_int_n))
        }
        
        dat_cox <- create_cox_dat(path_list = path_list, t_star = t_star)
        dat_cox <- dat_cox %>% left_join(df_var, by = "ID")
        
        col_ind <- which(colnames(dat_cox) %in% c("Treat", "Relapse", "Var1", "Var2"))
        dat_cox[, col_ind] <- lapply(dat_cox[, col_ind], function(x) factor(x, levels = c(0, 1)))
        
        # Avoid zero-length interval between start and stop due to high precision data      
        threshold <- 1e-06
        for (i in seq_len(nrow(dat_cox))) {
          if (abs(dat_cox$Start[i] - dat_cox$Stop[i]) < threshold) {
            id_data <- dat_cox[dat_cox$ID == dat_cox$ID[i], ]
            id_index <- which(rownames(id_data) == rownames(dat_cox[i, ]))
            
            if (id_index == nrow(id_data)) {
              dat_cox$Stop[i] <- dat_cox$Start[i] + threshold
            } else {
              dat_cox$Stop[i] <- dat_cox$Start[i] + threshold
              dat_cox$Start[i + 1] <- dat_cox$Start[i] + threshold
            }
          }
        }
        
        formula <- as.formula(cox_fit_list[[model]])
        cox_fit <- coxph(formula, data = dat_cox)
        
        # 1. model coefficients / discrete hazard ratio
        cox_coef <- summary(cox_fit)$coefficients
        
        # 2. ratio of cumulative hazard
        # 2.1 baseline cumulative hazard function
        if (scale == "forward" & model == 1) {
          baseline_dat <- data.frame(Treat = factor(0, levels = c(0, 1)), 
                                     Relapse = factor(0, levels = c(0, 1)), 
                                     Var1 = factor(0, levels = c(0, 1)), 
                                     Var2 = factor(0, levels = c(0, 1)))
        } else {
          baseline_dat <- data.frame(Treat = factor(0, levels = c(0, 1)), 
                                     Relapse = factor(0, levels = c(0, 1)), 
                                     Var1 = factor(0, levels = c(0, 1)), 
                                     Var2 = factor(0, levels = c(0, 1)), 
                                     Start = 0)
        }
        
        base_haz <- basehaz(cox_fit, newdata = baseline_dat)
        base_haz <- rbind(c(0, 0), base_haz)
        haz_fun <- splinefun(base_haz$time, base_haz$hazard, method = "monoH.FC")
        
        # 2.2 calculate ratio of cumulative hazards, in reference to never treat
        # 2.2.1 never treat (reference group)
        dat_cox_int_0$exp_lp <- predict(cox_fit, type = "risk", newdata = dat_cox_int_0)
        
        # calculate step-wise cumulative hazard
        dat_cox_int_0$discret_haz <- haz_fun(dat_cox_int_0$Stop) - haz_fun(dat_cox_int_0$Start)
        dat_cox_int_0$cumhaz_step <- dat_cox_int_0$discret_haz * dat_cox_int_0$exp_lp
        
        cumhaz_0 <- unlist(tapply(dat_cox_int_0$cumhaz_step, dat_cox_int_0$ID, FUN = sum))
        
        est_abs_risk <- mean(1 - exp(-cumhaz_0))
        
        # 2.2.2 treat under specific intervention strategies
        est_cumhaz_contrast <- numeric()
        est_risk_contrast <- numeric()
        
        for (int_time in int_list) {
          dat_cox_int_1 <- dat_cox_int_n_list[[which(int_list == int_time)]]
          dat_cox_int_1$exp_lp <- predict(cox_fit, type = "risk", newdata = dat_cox_int_1)
          
          # calculate step-wise cumulative hazard
          dat_cox_int_1$discret_haz <- haz_fun(dat_cox_int_1$Stop) - haz_fun(dat_cox_int_1$Start)
          dat_cox_int_1$cumhaz_step <- dat_cox_int_1$discret_haz * dat_cox_int_1$exp_lp
          
          cumhaz_1 <- unlist(tapply(dat_cox_int_1$cumhaz_step, dat_cox_int_1$ID, FUN = sum))
          
          est_abs_risk <- c(est_abs_risk, mean(1 - exp(-cumhaz_1)))
          est_cumhaz_contrast[which(int_list == int_time)] <- mean(cumhaz_1)/mean(cumhaz_0)
          est_risk_contrast[which(int_list == int_time)] <- mean(1 - exp(-cumhaz_1))/mean(1 - exp(-cumhaz_0))
        }
        
        cox_results <- list(cox_coef = cox_coef, 
                            est_abs_risk = est_abs_risk,
                            est_cumhaz_contrast = est_cumhaz_contrast,
                            est_risk_contrast = est_risk_contrast)
        
        return(cox_results)
      }
      
      if (scale == "forward") {
        file_name <- paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_fit_batch_", batch, ".rds")
      } else {
        file_name <- switch(model,
                            "1" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_linear_batch_", batch, ".rds"),
                            "2" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_quadratic_batch_", batch, ".rds"),
                            "3" = paste0("./model_fit/", scale, "/cox_td/", ms_scale, "_cox_pspline_batch_", batch, ".rds"))
      }
      
      saveRDS(cox_results_batch, file = file_name)
    }
  }
}

stopCluster(cl)
