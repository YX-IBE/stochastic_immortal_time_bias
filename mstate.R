# Multistate modelling

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
            paste0("./model_fit/", scale, "/mstate"))
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir)
      print(paste("Directory", dir, "created."))
    } else {
      print(paste("Directory", dir, "already exists."))
    }
  }
  
  # Load each batch sequentially
  for (batch in 1 : n_batch) {
    
    file_name <- paste0("./model_fit/", scale, "/mstate/", ms_scale, "_batch_", batch, ".rds")
    
    # Skip processing if the file already exists
    if (file.exists(file_name)) {
      cat(sprintf("File already exists: %s. Skipping...\n", file_name))
      next
    }
    
    path_lists_batch <- readRDS(paste0("./simulation/", scale, "/path_lists_batch_", batch, ".rds"))
    
    # fit the model in parallel for each path_list in the current batch
    mstate_results_batch <- foreach(simul = 1 : batch_size, .combine = function(x, y) c(x, list(y)), .init = list(), .packages = c("tidyverse", "survival", "mstate")) %dopar% {
      source("functions.R")
      
      path_list <- path_lists_batch[[simul]]
      
      # fit multi-state model
      dat_mstate_prep <- create_wide_mstate_dat(path_list = path_list, t_star = t_star) 
      dat_mstate_prep <- dat_mstate_prep %>% left_join(df_var, by = "ID")
      
      col_ind <- which(colnames(dat_mstate_prep) %in% c("Var1", "Var2"))
      dat_mstate_prep[, col_ind] <- lapply(dat_mstate_prep[, col_ind], function(x) factor(x, levels = c(0, 1)))
      
      tmat <- matrix(NA, 8, 8)
      tmat[1, c(2, 3, 5)] <- 1:3
      tmat[2, 4] <- 4
      tmat[5, 6:7] <- 5:6
      tmat[6, 8] <- 7
      dimnames(tmat) <- list(from = c("000", "010", "001", "011", "100", "110", "101", "111"), 
                             to = c("000", "010", "001", "011", "100", "110", "101", "111"))
      
      covs_mstate <- c("Var1", "Var2", "t_2", "t_5", "t_6")
      
      dat_mstate <- msprep(time = c(NA, "t_2", "t_3", "t_4", "t_5", "t_6", "t_7", "t_8"), 
                           status = c(NA, "s_2", "s_3", "s_4", "s_5", "s_6", "s_7", "s_8"), 
                           data = dat_mstate_prep, trans = tmat, keep = covs_mstate)
      dat_mstate <- expand.covs(dat_mstate, covs_mstate, append = T, longnames = F)
      
      # Avoid zero-length interval between start and stop due to high precision data      
      threshold <- 1e-06
      for (i in seq_len(nrow(dat_mstate))) {
        if (abs(dat_mstate$Tstart[i] - dat_mstate$Tstop[i]) < threshold) {
          id_data <- dat_mstate[dat_mstate$id == dat_mstate$id[i], ]
          index <- which(id_data$time < threshold)
          
          id_data$Tstop[index] <- id_data$Tstop[index] + threshold
          id_data$time[index] <- id_data$time[index] + threshold
          
          if (any(id_data$from %in% id_data$to[index])) {
            id_data$Tstart[id_data$from %in% id_data$to[index]] <- id_data$Tstart[id_data$from %in% id_data$to[index]] + threshold
          }
          dat_mstate[dat_mstate$id == dat_mstate$id[i], ] <- id_data
        }
      }
      
      dat_mstate$strata <- ifelse(dat_mstate$trans == 3, 1,
                                  ifelse(dat_mstate$trans %in% c(1, 5), 2, 3))
      
      dat_mstate$his_death <- 0
      dat_mstate$his_death[dat_mstate$trans == 4] <- 1
      dat_mstate$his_death[dat_mstate$trans == 6] <- 2
      dat_mstate$his_death[dat_mstate$trans == 7] <- 3
      dat_mstate$his_death <- factor(dat_mstate$his_death)
      
      dat_mstate$his_rel <- 0
      dat_mstate$his_rel[dat_mstate$trans == 5] <- 1
      dat_mstate$his_rel <- factor(dat_mstate$his_rel)
      
      if (ms_scale == "forward") {
        mstate_fit <- coxph(Surv(Tstart, Tstop, status) ~ 
                              Var1.1 + Var1.2 + Var1.3 + Var1.4 + Var1.5 + Var1.6 + Var1.7 + 
                              Var2.1 + Var2.2 + Var2.3 + Var2.4 + Var2.5 + Var2.6 + Var2.7 + 
                              t_2.4 + I(t_2.4^2) + t_5.5 + I(t_5.5^2) + 
                              t_5.6 + I(t_5.6^2) + t_6.7 + I(t_6.7^2) + his_death + his_rel +
                              strata(strata), data = dat_mstate, method = "breslow")
      } else {
        mstate_fit <- coxph(Surv(time, status) ~ 
                              Var1.1 + Var1.2 + Var1.3 + Var1.4 + Var1.5 + Var1.6 + Var1.7 + 
                              Var2.1 + Var2.2 + Var2.3 + Var2.4 + Var2.5 + Var2.6 + Var2.7 + 
                              t_2.4 + I(t_2.4^2) + t_5.5 + I(t_5.5^2) + 
                              t_5.6 + I(t_5.6^2) + t_6.7 + I(t_6.7^2) + his_death + his_rel +
                              strata(strata), data = dat_mstate, method = "breslow")
      }
      
      # 1. model coefficients
      mstate_coef <- summary(mstate_fit)$coefficients
      
      # 2.1 cumulative hazard if never treat (reference group)
      
      sim_dat_int_0 <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_36_", ms_scale, ".rds"))
      sim_dat_int_0$Treat_s <- 0
      
      sim_dat_int_0$s_2 <- sim_dat_int_0$Relapse_s
      sim_dat_int_0$t_2 <- sim_dat_int_0$Relapse_s * sim_dat_int_0$Relapse_t
      sim_dat_int_0$s_5 <- 0
      sim_dat_int_0$t_5 <- 0
      sim_dat_int_0$s_6 <- 0
      sim_dat_int_0$t_6 <- 0
      
      cumhaz_0 <- NULL
      
      for (i in seq_len(nrow(sim_dat_int_0))) {
        df_i <- sim_dat_int_0[i, ]
        
        newd_i <- data.frame(Var1 = rep(df_i$Var1, 7), Var2 = rep(df_i$Var2, 7), 
                             t_2 = rep(df_i$t_2, 7), 
                             t_5 = rep(df_i$t_5, 7), 
                             t_6 = rep(df_i$t_6, 7), 
                             trans = 1:7)
        
        attr(newd_i, "trans") <- tmat
        class(newd_i) <- c("msdata", "data.frame")
        newd_i <- expand.covs(newd_i, covs_mstate, longnames = F)
        
        newd_i$strata <- ifelse(newd_i$trans == 3, 1,
                                ifelse(newd_i$trans %in% c(1, 5), 2, 3))
        
        newd_i$his_death <- 0
        newd_i$his_death[newd_i$trans == 4] <- 1
        newd_i$his_death[newd_i$trans == 6] <- 2
        newd_i$his_death[newd_i$trans == 7] <- 3
        newd_i$his_death <- factor(newd_i$his_death)
        
        newd_i$his_rel <- 0
        newd_i$his_rel[newd_i$trans == 5] <- 1
        newd_i$his_rel <- factor(newd_i$his_rel)
        
        msfit_haz <- msfit(object = mstate_fit, newdata = newd_i, trans = tmat)
        
        haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 2, ]
        haz_dat <- rbind(c(0, 0, 2), haz_dat)
        haz1 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
        
        if (df_i$s_2 == 1) {
          haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 4, ]
          haz_dat <- rbind(c(0, 0, 4), haz_dat)
          haz2 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
          
          cumhaz_0_i <- haz2(t_star) - haz2(df_i$t_2) + haz1(df_i$t_2)
          
        } else {
          cumhaz_0_i <- haz1(t_star)
        }
        
        cumhaz_0 <- c(cumhaz_0, cumhaz_0_i)
      }
      
      est_abs_risk <- mean(1 - exp(-cumhaz_0))
      
      # 2.2 cumulative hazard if under specific intervention strategies
      est_cumhaz_contrast <- numeric()
      est_risk_contrast <- numeric()
      
      for (int_time in int_list) {
        
        sim_dat_int_1 <- readRDS(paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_", int_time, "_", ms_scale, ".rds"))
        
        sim_dat_int_1$s_2 <- as.numeric(sim_dat_int_1$Relapse_t < sim_dat_int_1$Treat_t)
        sim_dat_int_1$t_2 <- sim_dat_int_1$s_2 * sim_dat_int_1$Relapse_t
        
        sim_dat_int_1$s_5 <- sim_dat_int_1$Treat_s
        sim_dat_int_1$t_5 <- sim_dat_int_1$Treat_s * sim_dat_int_1$Treat_t
        
        sim_dat_int_1$s_6 <- as.numeric(sim_dat_int_1$Relapse_t >= sim_dat_int_1$Treat_t & sim_dat_int_1$Relapse_t < t_star)
        sim_dat_int_1$t_6 <- sim_dat_int_1$s_6 * sim_dat_int_1$Relapse_t
        
        cumhaz_1 <- NULL
        
        for (i in seq_len(nrow(sim_dat_int_1))) {
          df_i <- sim_dat_int_1[i, ]
          
          newd_i <- data.frame(Var1 = rep(df_i$Var1, 7), Var2 = rep(df_i$Var2, 7), 
                               t_2 = rep(df_i$t_2, 7), 
                               t_5 = rep(df_i$t_5, 7), 
                               t_6 = rep(df_i$t_6, 7), 
                               trans = 1:7)
          
          attr(newd_i, "trans") <- tmat
          class(newd_i) <- c("msdata", "data.frame")
          newd_i <- expand.covs(newd_i, covs_mstate, longnames = F)
          
          newd_i$strata <- ifelse(newd_i$trans == 3, 1,
                                  ifelse(newd_i$trans %in% c(1, 5), 2, 3))
          
          newd_i$his_death <- 0
          newd_i$his_death[newd_i$trans == 4] <- 1
          newd_i$his_death[newd_i$trans == 6] <- 2
          newd_i$his_death[newd_i$trans == 7] <- 3
          newd_i$his_death <- factor(newd_i$his_death)
          
          newd_i$his_rel <- 0
          newd_i$his_rel[newd_i$trans == 5] <- 1
          newd_i$his_rel <- factor(newd_i$his_rel)
          
          msfit_haz <- msfit(object = mstate_fit, newdata = newd_i, trans = tmat)
          
          haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 2, ]
          haz_dat <- rbind(c(0, 0, 2), haz_dat)
          haz1 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
          
          if (df_i$s_2 == 1) {
            haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 4, ]
            haz_dat <- rbind(c(0, 0, 4), haz_dat)
            haz2 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
            
            cumhaz_1_i <- haz2(t_star) - haz2(df_i$t_2) + haz1(df_i$t_2)
            
          } else if (df_i$s_6 == 0) {
            haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 6, ]
            haz_dat <- rbind(c(0, 0, 6), haz_dat)
            haz2 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
            
            cumhaz_1_i <- haz2(t_star) - haz2(df_i$t_5) + haz1(df_i$t_5)
            
          } else {
            haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 6, ]
            haz_dat <- rbind(c(0, 0, 6), haz_dat)
            haz2 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
            
            haz_dat <- msfit_haz$Haz[msfit_haz$Haz$trans == 7, ]
            haz_dat <- rbind(c(0, 0, 7), haz_dat)
            haz3 <- splinefun(haz_dat$time, haz_dat$Haz, method = "monoH.FC")
            
            cumhaz_1_i <- haz3(t_star) - haz3(df_i$t_6) + 
              haz2(df_i$t_6) - haz2(df_i$t_5) + 
              haz1(df_i$t_5)
          }
          
          cumhaz_1 <- c(cumhaz_1, cumhaz_1_i)
        }        
        
        est_abs_risk <- c(est_abs_risk, mean(1 - exp(-cumhaz_1)))
        est_cumhaz_contrast[which(int_list == int_time)] <- mean(cumhaz_1)/mean(cumhaz_0)
        est_risk_contrast[which(int_list == int_time)] <- mean(1 - exp(-cumhaz_1))/mean(1 - exp(-cumhaz_0))
      }
      
      mstate_results <- list(mstate_coef = mstate_coef, 
                             est_abs_risk = est_abs_risk,
                             est_cumhaz_contrast = est_cumhaz_contrast,
                             est_risk_contrast = est_risk_contrast)
      
      return(mstate_results)
    }
    
    file_name <- paste0("./model_fit/", scale, "/mstate/", ms_scale, "_batch_", batch, ".rds")
    saveRDS(mstate_results_batch, file = file_name)
  }
}

stopCluster(cl)
