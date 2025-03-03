setwd("~/S1")

library(foreach)
library(doParallel)
library(tidyverse)
library(survival)
library(mstate)

# Setup parallel backend
Sys.setenv(OPENBLAS_NUM_THREADS = 2)
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
int_list <- c(1, 3, 6, 9, 12, 36)
sim_pat <- 1000
RNG <- 1809

# Var1 and Var2
df_var <- as.data.frame(cbind(ID = 1:(n_pat * 4), Var1 = rep(NA, (n_pat * 4)), Var2 = rep(NA, (n_pat * 4))))
df_var$Var1 <- rep(rep(c(0, 1), 2), n_pat)
df_var$Var2 <- rep(rep(c(0, 1), c(2, 2)), n_pat)
df_var[, c("Var1", "Var2")] <- lapply(df_var[, c("Var1", "Var2")], function(x) factor(x, levels = c(0, 1)))

# we use mstate to estimate Aalen–Johansen-type estimator
# transition matrix
tmat <- matrix(NA, 8, 8)
tmat[1, c(2, 3, 5)] <- 1:3
tmat[2, 4] <- 4
tmat[5, 6:7] <- 5:6
tmat[6, 8] <- 7
dimnames(tmat) <- list(from = c("000", "010", "001", "011", "100", "110", "101", "111"), 
                       to = c("000", "010", "001", "011", "100", "110", "101", "111"))

mstate_covs <- c("Var1", "Var2", "t_2", "t_5", "t_6")

# Aalen-Johansen-type estimator is only for generating hypothetical population; where trans 1 and trans 5 are relevant
# trans 1: 000 to 010; 001 and 100 are competing events
# trans 5: 100 to 110; 101 are competing events
# only t_5 (time of treatment) will impact baseline hazard estimate for trans 5
newd_list <- list()
for (int_time in int_list) {
  newd_1 <- data.frame(Var1 = rep(0, 7), Var2 = rep(0, 7), 
                       t_2 = rep(0, 7), t_5 = rep(int_time, 7), t_6 = rep(0, 7), 
                       trans = 1:7)
  
  newd_1$Var1 <- factor(newd_1$Var1, levels = c(0, 1))
  newd_1$Var2 <- factor(newd_1$Var2, levels = c(0, 1))
  
  attr(newd_1, "trans") <- tmat
  class(newd_1) <- c("msdata", "data.frame")
  newd_1 <- expand.covs(newd_1, mstate_covs, longnames = F)
  
  newd_1$strata <- ifelse(newd_1$trans == 3, 1,
                          ifelse(newd_1$trans %in% c(1, 5), 2, 3))
  
  newd_1$his_death <- 0
  newd_1$his_death[newd_1$trans == 4] <- 1
  newd_1$his_death[newd_1$trans == 6] <- 2
  newd_1$his_death[newd_1$trans == 7] <- 3
  newd_1$his_death <- factor(newd_1$his_death)
  
  newd_1$his_rel <- 0
  newd_1$his_rel[newd_1$trans == 5] <- 1
  newd_1$his_rel <- factor(newd_1$his_rel)
  
  newd_2 <- data.frame(Var1 = rep(1, 7), Var2 = rep(0, 7), 
                       t_2 = rep(0, 7), t_5 = rep(int_time, 7), t_6 = rep(0, 7), 
                       trans = 1:7)
  
  newd_2$Var1 <- factor(newd_2$Var1, levels = c(0, 1))
  newd_2$Var2 <- factor(newd_2$Var2, levels = c(0, 1))
  
  attr(newd_2, "trans") <- tmat
  class(newd_2) <- c("msdata", "data.frame")
  newd_2 <- expand.covs(newd_2, mstate_covs, longnames = F)
  
  newd_2$strata <- ifelse(newd_2$trans == 3, 1,
                          ifelse(newd_2$trans %in% c(1, 5), 2, 3))
  
  newd_2$his_death <- 0
  newd_2$his_death[newd_2$trans == 4] <- 1
  newd_2$his_death[newd_2$trans == 6] <- 2
  newd_2$his_death[newd_2$trans == 7] <- 3
  newd_2$his_death <- factor(newd_2$his_death)
  
  newd_2$his_rel <- 0
  newd_2$his_rel[newd_2$trans == 5] <- 1
  newd_2$his_rel <- factor(newd_2$his_rel)
  
  newd_3 <- data.frame(Var1 = rep(0, 7), Var2 = rep(1, 7), 
                       t_2 = rep(0, 7), t_5 = rep(int_time, 7), t_6 = rep(0, 7), 
                       trans = 1:7)
  
  newd_3$Var1 <- factor(newd_3$Var1, levels = c(0, 1))
  newd_3$Var2 <- factor(newd_3$Var2, levels = c(0, 1))
  
  attr(newd_3, "trans") <- tmat
  class(newd_3) <- c("msdata", "data.frame")
  newd_3 <- expand.covs(newd_3, mstate_covs, longnames = F)
  
  newd_3$strata <- ifelse(newd_3$trans == 3, 1,
                          ifelse(newd_3$trans %in% c(1, 5), 2, 3))
  
  newd_3$his_death <- 0
  newd_3$his_death[newd_3$trans == 4] <- 1
  newd_3$his_death[newd_3$trans == 6] <- 2
  newd_3$his_death[newd_3$trans == 7] <- 3
  newd_3$his_death <- factor(newd_3$his_death)
  
  newd_3$his_rel <- 0
  newd_3$his_rel[newd_3$trans == 5] <- 1
  newd_3$his_rel <- factor(newd_3$his_rel)
  
  newd_4 <- data.frame(Var1 = rep(1, 7), Var2 = rep(1, 7), 
                       t_2 = rep(0, 7), t_5 = rep(int_time, 7), t_6 = rep(0, 7), 
                       trans = 1:7)
  
  newd_4$Var1 <- factor(newd_4$Var1, levels = c(0, 1))
  newd_4$Var2 <- factor(newd_4$Var2, levels = c(0, 1))
  
  attr(newd_4, "trans") <- tmat
  class(newd_4) <- c("msdata", "data.frame")
  newd_4 <- expand.covs(newd_4, mstate_covs, longnames = F)
  
  newd_4$strata <- ifelse(newd_4$trans == 3, 1,
                          ifelse(newd_4$trans %in% c(1, 5), 2, 3))
  
  newd_4$his_death <- 0
  newd_4$his_death[newd_4$trans == 4] <- 1
  newd_4$his_death[newd_4$trans == 6] <- 2
  newd_4$his_death[newd_4$trans == 7] <- 3
  newd_4$his_death <- factor(newd_4$his_death)
  
  newd_4$his_rel <- 0
  newd_4$his_rel[newd_4$trans == 5] <- 1
  newd_4$his_rel <- factor(newd_4$his_rel)
  
  newd_list[[which(int_list == int_time)]] <- list(newd_1, newd_2, newd_3, newd_4)
}

for (scale in scales) {
  
  # create a folder for each scale
  dirs <- c("./history_int",
            paste0("./history_int/", scale))
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      print(paste("Directory", dir, "created."))
    } else {
      print(paste("Directory", dir, "already exists."))
    }
  }
  
  # loop over 10 batch
  # the 100 tasks within each batch is computed in parallel
  for (batch in 1 : n_batch) {
    path_lists_batch <- readRDS(paste0("./simulation/", scale, "/path_lists_batch_", batch, ".rds"))
    batch_rng <- RNG + batch * 100
    
    foreach(simul = 1 : batch_size, .packages = c("tidyverse", "mstate"), .options.RNG = batch_rng) %dopar% {
      source("functions.R")
      path_list <- path_lists_batch[[simul]]
      
      # prepare data for multistate model
      mstate_dat_wide <- create_wide_mstate_dat(path_list, t_star)
      mstate_dat_wide <- mstate_dat_wide %>% left_join(df_var, by = "ID")
      
      mstate_covs <- c("Var1", "Var2", "t_2", "t_5", "t_6")
      mstate_dat <- msprep(time = c(NA, "t_2", "t_3", "t_4", "t_5", "t_6", "t_7", "t_8"), 
                           status = c(NA, "s_2", "s_3", "s_4", "s_5", "s_6", "s_7", "s_8"), 
                           data = mstate_dat_wide, trans = tmat, keep = mstate_covs)
      mstate_dat <- expand.covs(mstate_dat, mstate_covs, append = T, longnames = F)
      
      # Avoid zero-length interval between start and stop due to high precision data      
      threshold <- 1e-06
      for (i in seq_len(nrow(mstate_dat))) {
        if (abs(mstate_dat$Tstart[i] - mstate_dat$Tstop[i]) < threshold) {
          id_data <- mstate_dat[mstate_dat$id == mstate_dat$id[i], ]
          index <- which(id_data$time < threshold)
          
          id_data$Tstop[index] <- id_data$Tstop[index] + threshold
          id_data$time[index] <- id_data$time[index] + threshold
          
          if (any(id_data$from %in% id_data$to[index])) {
            id_data$Tstart[id_data$from %in% id_data$to[index]] <- id_data$Tstart[id_data$from %in% id_data$to[index]] + threshold
          }
          mstate_dat[mstate_dat$id == mstate_dat$id[i], ] <- id_data
        }
      }
      
      mstate_dat$strata <- ifelse(mstate_dat$trans == 3, 1,
                                  ifelse(mstate_dat$trans %in% c(1, 5), 2, 3))
      
      mstate_dat$his_death <- 0
      mstate_dat$his_death[mstate_dat$trans == 4] <- 1
      mstate_dat$his_death[mstate_dat$trans == 6] <- 2
      mstate_dat$his_death[mstate_dat$trans == 7] <- 3
      mstate_dat$his_death <- factor(mstate_dat$his_death)
      
      mstate_dat$his_rel <- 0
      mstate_dat$his_rel[mstate_dat$trans == 5] <- 1
      mstate_dat$his_rel <- factor(mstate_dat$his_rel)
      
      # fit with two scales
      mstate_fit_forward <- coxph(Surv(Tstart, Tstop, status) ~ 
                                    Var1.1 + Var1.2 + Var1.3 + Var1.4 + Var1.5 + Var1.6 + Var1.7 + 
                                    Var2.1 + Var2.2 + Var2.3 + Var2.4 + Var2.5 + Var2.6 + Var2.7 + 
                                    t_2.4 + I(t_2.4^2) + t_5.5 + I(t_5.5^2) + 
                                    t_5.6 + I(t_5.6^2) + t_6.7 + I(t_6.7^2) + his_death + his_rel +
                                    strata(strata), data = mstate_dat, method = "breslow")
      
      mstate_fit_renewal <- coxph(Surv(time, status) ~ 
                                    Var1.1 + Var1.2 + Var1.3 + Var1.4 + Var1.5 + Var1.6 + Var1.7 + 
                                    Var2.1 + Var2.2 + Var2.3 + Var2.4 + Var2.5 + Var2.6 + Var2.7 + 
                                    t_2.4 + I(t_2.4^2) + t_5.5 + I(t_5.5^2) + 
                                    t_5.6 + I(t_5.6^2) + t_6.7 + I(t_6.7^2) + his_death + his_rel +
                                    strata(strata), data = mstate_dat, method = "breslow")
      
      ms_scale <- ifelse(scale == "extended", "renewal", scale)
      
      # loop over intervention strategies for the four subgroups, each with 1000 participants
      for (int_time in int_list) {
        
        newd_list_int <- newd_list[[which(int_list == int_time)]]
        
        dat_int_sim <- NULL
        
        for (group in 1:4) {
          
          # fit the Aalen–Johansen type estimator for each subgroup
          msfit_haz <- msfit(object = get(paste0("mstate_fit_", ms_scale)), newdata = newd_list_int[[group]], trans = tmat)
          
          # cumulative hazard for relapse from 100 to 110 over (0, t_star)
          trt1 <- msfit_haz$Haz[msfit_haz$Haz$trans == 5, ]
          trt1 <- rbind(c(0, 0, 5), trt1)
          if (abs(trt1$time[nrow(trt1)] - t_star) >= 1e-6) {
            trt1 <- rbind(trt1, c(t_star, trt1$Haz[nrow(trt1)], 5))
          }
          haz_rel_1 <- splinefun(trt1$time, trt1$Haz, method = "monoH.FC")
          
          # cumulative hazard for relapse from 000 to 010 over (0, t_star)
          trt0 <- msfit_haz$Haz[msfit_haz$Haz$trans == 1, ]
          trt0 <- rbind(c(0, 0, 1), trt0)
          if (abs(trt0$time[nrow(trt0)] - t_star) >= 1e-6) {
            trt0 <- rbind(trt0, c(t_star, trt0$Haz[nrow(trt0)], 1))
          }
          haz_rel_0 <- splinefun(trt0$time, trt0$Haz, method = "monoH.FC")
          
          if (int_time < t_star) {
            
            # CDF for relapse if treated at int_time
            cdf_1 <- function(t) {
              ifelse(t <= int_time,
                     1 - exp(-haz_rel_0(t)), 
                     1 - exp(-(haz_rel_1(t) - haz_rel_1(int_time) + haz_rel_0(int_time)))
              )
            }
            
            t_seq <- seq(0, t_star, by = 0.1)
            cdf_1_values <- sapply(t_seq, cdf_1)
            
            # adjust if more than 3 consecutive cdf_1_values stabilize
            consecutive_n <- 4
            diff_values <- diff(cdf_1_values)
            stable_point <- NA
            
            # check from the second time point
            for (i in 2 : (length(diff_values) - consecutive_n + 1)) {
              if (all(abs(diff_values[i : (i + consecutive_n - 1)]) < 1e-6)) {
                stable_point <- i
                break
              }
            }
            
            if (!is.na(stable_point)) {
              cdf_1_values_adjusted <- cdf_1_values
              stable_t <- t_seq[stable_point + 1]
              adjusted_t_seq <- seq(stable_t, t_star, by = 0.1)
              adjusted_values <- 1e-6 * seq_along(adjusted_t_seq)
              cdf_1_values_adjusted[(stable_point + 1) : length(cdf_1_values_adjusted)] <- cdf_1_values_adjusted[(stable_point + 1) : length(cdf_1_values_adjusted)] + adjusted_values
              cdf_1_values <- cdf_1_values_adjusted
            }
            
            spar <- 0.3
            by_int <- 0.1
            success <- FALSE
            
            while (!success) {
              tryCatch({
                smooth_fit <- smooth.spline(t_seq, cdf_1_values, spar = spar)
                smooth_t_seq <- seq(0, t_star, by = by_int)
                smoothed_cdf_1_values <- predict(smooth_fit, smooth_t_seq)$y
                smoothed_cdf_1_values <- pmax(cummax(smoothed_cdf_1_values), 0)
                smoothed_cdf_1_values[1] <- 0
                smoothed_cdf_1_values[length(smoothed_cdf_1_values)] <- max(cdf_1_values[length(cdf_1_values)], smoothed_cdf_1_values[length(smoothed_cdf_1_values)])
                cdf_1_smoothed <- splinefun(smooth_t_seq, smoothed_cdf_1_values, method = "monoH.FC")
                
                lambda <- 0.05
                cdf_1_extended <- function(t) {
                  ifelse(t <= t_star,
                         cdf_1_smoothed(t),
                         cdf_1_smoothed(t_star) + (1 - cdf_1_smoothed(t_star)) * (1 - exp(-lambda * (t - t_star)))
                  )
                }
                
                gen <- pinv.new(cdf = cdf_1_extended, lb = 0, ub = Inf, uresolution = 1e-12)
                success <- TRUE
              }, error = function(e) {
                if (by_int >= 2) {
                  spar <<- spar + 0.1
                } else {
                  by_int <<- by_int + 0.1
                }
                
                if (spar >= 1) {
                  write(NA, file = paste0("./history_int/", "failed_file_", scale, "_batch_", batch, "_simul_", simul, "_int_", int_time,  "_group_", group, "_treat.txt"), append = TRUE)
                  success <<- TRUE
                }
              })
            }
            
            if (success) {
              dat_int_sim <- rbind(dat_int_sim, cbind(uq(gen, runif(sim_pat, 0, 1)), group, spar, by_int, Treat_s = 1, Treat_t = int_time))
            }
            
          } else {

            # CDF for relapse if never treat
            cdf_0 <- function(t) {
              1 - exp(-haz_rel_0(t))
            }
            
            t_seq <- seq(0, t_star, by = 0.1)
            cdf_0_values <- sapply(t_seq, cdf_0)
            
            # adjust if more than 3 consecutive cdf_0_values stabilize
            consecutive_n <- 4
            diff_values <- diff(cdf_0_values)
            stable_point <- NA
            
            # check from the second time point
            for (i in 2 : (length(diff_values) - consecutive_n + 1)) {
              if (all(abs(diff_values[i : (i + consecutive_n - 1)]) < 1e-6)) {
                stable_point <- i
                break
              }
            }
            
            if (!is.na(stable_point)) {
              cdf_0_values_adjusted <- cdf_0_values
              stable_t <- t_seq[stable_point + 1]
              adjusted_t_seq <- seq(stable_t, t_star, by = 0.1)
              adjusted_values <- 1e-6 * seq_along(adjusted_t_seq)
              cdf_0_values_adjusted[(stable_point + 1) : length(cdf_0_values_adjusted)] <- cdf_0_values_adjusted[(stable_point + 1) : length(cdf_0_values_adjusted)] + adjusted_values
              cdf_0_values <- cdf_0_values_adjusted
            }
            
            spar <- 0.3
            by_int <- 0.1
            success <- FALSE
            
            while (!success) {
              tryCatch({
                smooth_fit <- smooth.spline(t_seq, cdf_0_values, spar = spar)
                smooth_t_seq <- seq(0, t_star, by = by_int)
                smoothed_cdf_0_values <- predict(smooth_fit, smooth_t_seq)$y
                smoothed_cdf_0_values <- pmax(cummax(smoothed_cdf_0_values), 0)
                smoothed_cdf_0_values[1] <- 0
                smoothed_cdf_0_values[length(smoothed_cdf_0_values)] <- max(cdf_0_values[length(cdf_0_values)], smoothed_cdf_0_values[length(smoothed_cdf_0_values)])
                cdf_0_smoothed <- splinefun(smooth_t_seq, smoothed_cdf_0_values, method = "monoH.FC")
                
                lambda <- 0.05
                cdf_0_extended <- function(t) {
                  ifelse(t <= t_star,
                         cdf_0_smoothed(t),
                         cdf_0_smoothed(t_star) + (1 - cdf_0_smoothed(t_star)) * (1 - exp(-lambda * (t - t_star)))
                  )
                }
                
                gen <- pinv.new(cdf = cdf_0_extended, lb = 0, ub = Inf, uresolution = 1e-12)
                success <- TRUE
              }, error = function(e) {
                if (by_int >= 2) {
                  spar <<- spar + 0.1
                } else {
                  by_int <<- by_int + 0.1
                }
                
                if (spar >= 1) {
                  write(NA, file = paste0("./history_int/", "failed_file_", scale, "_batch_", batch, "_simul_", simul, "_int_", int_time, "_group_", group, "_notreat.txt"), append = TRUE)
                  success <<- TRUE
                }
              })
            }
            
            if (success) {
              dat_int_sim <- rbind(dat_int_sim, cbind(uq(gen, runif(sim_pat, 0, 1)), group, spar, by_int, Treat_s = 0, Treat_t = int_time))
            }
            
          }
        }
        dat_int_sim <- as.data.frame(dat_int_sim)
        dat_int_sim$Var1 <- factor(ifelse(dat_int_sim$group %in% c(1, 3), 0, 1))
        dat_int_sim$Var2 <- factor(ifelse(dat_int_sim$group %in% c(1, 2), 0, 1))
        dat_int_sim <- dat_int_sim[, -2]
        names(dat_int_sim)[1] <- "Relapse_t"
        dat_int_sim$Relapse_s <- ifelse(dat_int_sim$Relapse_t < t_star, 1, 0)
        if (int_time != t_star) {
          dat_int_sim$Treat_s <- ifelse(dat_int_sim$Relapse_t < dat_int_sim$Treat_t, 0, 1)
        }
        saveRDS(dat_int_sim, paste0("./history_int/", scale, "/batch_", batch, "_simul_", simul, "_int_", int_time, "_", ms_scale, ".rds"))
      }
    }
  }
}

stopCluster(cl)
