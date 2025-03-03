# Continuous Time Markov Jump processes
# Functions for Simulation
# Date: Fri 18.10.2024

# We simulate a continuous time Markov jump process, defined by three binary variables: treatment, relapse, and death
# Assume that: 
# (1) lambda_t and gamma_t relates to the time when treatment is administered
#     lambda_r and gamma_r relates to the time when relapse occurs, given no treatment
#     lambda_d and gamma_d relates to the time when death occurs, given no treatment or relapse
# (2) Treatment modifies time to relapse by hr_t_r
#     Treatment modifies time to death by hr_t_d
#     Relapse modifies time to death by hr_r_d
#     Effect modification is multiplicative, i.e. hr_t_d * hr_r_d
# (3) Treatment can be administered only before the occurrence of relapse and death

library(Runuran)
library(tidyverse)
library(dynpred)

# Assign matrix values
# matrix_vals[[1]]: hazard of the given distribution at time t, h(t)
# matrix_vals[[2]]: survival of the given distribution at time t, S(t)
# dist: character string specifying the distribution, built-in options: "exponential", "weibull", "gompertz", "llogistic"
assign_matrix_values <- function(t, dist, lambda, gamma = NA) {
  
  if (!dist %in% c("exponential", "weibull", "gompertz", "llogistic")) {
    stop("Current built-in distributions: exponential, Weibull, Gompertz, log-logistic")
  }
  
  if (!is.numeric(t) || !is.numeric(lambda)) {
    stop("t and lambda should be numeric values.")
  }
  
  if (dist != "exponential" && (is.na(gamma) || !is.numeric(gamma))) {
    stop("gamma must be provided as a numeric value for the specified distribution.")
  }
  
  if (dist == "exponential") {
    rate <- lambda
    surv <- exp(-lambda * t)
    
  } else if (dist == "weibull") {
    rate <- lambda * gamma * t^(gamma - 1)
    surv <- exp(-lambda * t^gamma)
    
  } else if (dist == "gompertz") {
    if (gamma == 0) {
      stop("Scale parameter gamma should not be 0 for Gompertz distribution: consider changing to exponential distribution.")
    }
    rate <- lambda * exp(gamma * t)
    surv <- exp(-lambda * (exp(gamma * t) - 1) / gamma)
    
  } else if (dist == "llogistic") {
    rate <- (lambda * gamma * t^(gamma - 1)) * (1 + lambda * t ^ gamma)^(-1)
    surv <- (1 + lambda * t^gamma) ^ (-1)
    
  } 
  
  matrix_vals <- list(rate = rate, surv = surv)
  return(matrix_vals)
}


# Create a list with two matrices
# mat_list[[1]]: generator matrix
# mat_list[[2]]: off-diagonal elements are survival probabilities at t; diagonal elements are set to 1
# dists = c(dist_t, dist_r, dist_d): a character vector specifying the parametric survival distributions
# lambdas = c(lambda_t, lambda_r, lambda_d): a numeric vector corresponding to the scale parameters
# gammas = c(gamma_t, gamma_r, gamma_d): a numeric vector corresponding to the shape parameters
# hr_t_r: hazard ratio of treatment on relapse
# hr_t_d: hazard ratio of treatment on death
# hr_r_d: hazard ratio of relapse on death
# relapse_his: TRUE/FALSE, does history prior to relapse, e.g. duration of complete remission, impact post-relapse hazard?
# his_var: a continuous variable used to account for history prior to relapse, e.g. duration of complete remission
# post_rel: TRUE/FALSE, should post relapse survival follow a new distribution?
# post_rel_dist_d: a character var specifying the new post-relapse parametric survival distribution
# post_rel_lambda_d: a numeric var corresponding to the scale parameter
# post_rel_gamma_d: a numeric var corresponding to the shape parameter
mat_list <- function(t, dists, lambdas, gammas,
                     hr_t_r, hr_t_d, hr_r_d,
                     relapse_his = F, his_var = NA,
                     post_rel = F, post_rel_dist_d = NA, 
                     post_rel_lambda_d = NA, post_rel_gamma_d = NA) {
  
  # Check if dists has exactly three elements for treatment, relapse, and death
  if (length(dists) != 3) {
    stop("dists should be specified individually as character vector for treatment, relapse, and death.")
  }
  
  if (!all(dists %in% c("exponential", "weibull", "gompertz", "llogistic")) || 
      (post_rel && !is.na(post_rel_dist_d) && !post_rel_dist_d %in% c("exponential", "weibull", "gompertz", "llogistic"))) {
    stop("Current built-in distributions: exponential, Weibull, Gompertz, log-logistic")
  }
  
  # Check if lambdas and gammas have the correct lengths
  if (length(lambdas) != 3 | length(gammas) != 3) {
    stop("lambdas and gammas should be specified individually as numeric vector for treatment, relapse, and death.
           Set gamma = NA or other arbitrary values as placeholder for exponential distribution.")
  }
  
  if (relapse_his & is.na(his_var)) {
    stop("relapse_his is TRUE, thus his_var should be specified.")
    
  } else if (relapse_his & is.numeric(hr_r_d)) {
    warning("hr_r_d is not time-dependent. To account for history prior to relapse under PH assumption (in relation to post-relapse baseline hazard), consider specifying a time-dependent hr_r_d as a string.")
    
  } else {
    hr_r_d_func <- function(his_var) eval(parse(text = hr_r_d))
  }
  
  if (post_rel) {
    if (is.na(post_rel_dist_d) || is.na(post_rel_lambda_d) || is.na(post_rel_gamma_d)) {
      stop("post_rel is TRUE, post_rel_dist_d, post_rel_lambda_d, and post_rel_gamma_d should all be specified.")
    }
  }
  
  g_mat <- matrix(0, ncol = 8, nrow = 8)
  s_mat <- matrix(1, ncol = 8, nrow = 8)
  colnames(g_mat) <- rownames(g_mat) <- colnames(s_mat) <- rownames(s_mat) <- c("000", "010", "001", "011", "100", "110", "101", "111")
  
  lambda_t <- lambdas[1]
  lambda_r <- lambdas[2]
  lambda_d <- lambdas[3]
  gamma_t <- gammas[1]
  gamma_r <- gammas[2]
  gamma_d <- gammas[3]
  dist_t <- dists[1]
  dist_r <- dists[2]
  dist_d <- dists[3]
  
  matrix_vals <- assign_matrix_values(t = t, dist = dist_r, lambda = lambda_r, gamma = gamma_r)
  g_mat[1, 2] <- matrix_vals$rate
  s_mat[1, 2] <- matrix_vals$surv
  
  matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d, gamma = gamma_d)
  g_mat[1, 3] <- matrix_vals$rate
  s_mat[1, 3] <- matrix_vals$surv
  
  matrix_vals <- assign_matrix_values(t = t, dist = dist_t, lambda = lambda_t, gamma = gamma_t)
  g_mat[1, 5] <- matrix_vals$rate
  s_mat[1, 5] <- matrix_vals$surv
  
  if (post_rel) {
    matrix_vals <- assign_matrix_values(t = t, dist = post_rel_dist_d, lambda = post_rel_lambda_d, gamma = post_rel_gamma_d)
    g_mat[2, 4] <- matrix_vals$rate
    s_mat[2, 4] <- matrix_vals$surv
    
  } else if (dist_d != "llogistic") {  
    if (exists("hr_r_d_func")) {
      matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d * hr_r_d_func(his_var), gamma = gamma_d)
    } else {
      matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d * hr_r_d, gamma = gamma_d)
    }
    g_mat[2, 4] <- matrix_vals$rate
    s_mat[2, 4] <- matrix_vals$surv
    
  } else {
    matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d, gamma = gamma_d)
    
    if (exists("hr_r_d_func")) {
      g_mat[2, 4] <- matrix_vals$rate * hr_r_d_func(his_var)
      s_mat[2, 4] <- exp(log(matrix_vals$surv) * hr_r_d_func(his_var))
      
    } else {
      g_mat[2, 4] <- matrix_vals$rate * hr_r_d
      s_mat[2, 4] <- exp(log(matrix_vals$surv) * hr_r_d) 
    }
  }
  
  matrix_vals <- assign_matrix_values(t = t, dist = dist_r, lambda = lambda_r * hr_t_r, gamma = gamma_r)
  g_mat[5, 6] <- matrix_vals$rate
  s_mat[5, 6] <- matrix_vals$surv
  
  if (dist_d != "llogistic") {  
    matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d * hr_t_d, gamma = gamma_d)
    g_mat[5, 7] <- matrix_vals$rate
    s_mat[5, 7] <- matrix_vals$surv
    
  } else {
    matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d, gamma = gamma_d)
    g_mat[5, 7] <- matrix_vals$rate * hr_t_d
    s_mat[5, 7] <- exp(log(matrix_vals$surv) * hr_t_d)
  }
  
  if (post_rel) {
    if (post_rel_dist_d != "llogistic") {
      matrix_vals <- assign_matrix_values(t = t, dist = post_rel_dist_d, lambda = post_rel_lambda_d * hr_t_d, gamma = post_rel_gamma_d)
      g_mat[6, 8] <- matrix_vals$rate
      s_mat[6, 8] <- matrix_vals$surv
      
    } else {
      matrix_vals <- assign_matrix_values(t = t, dist = post_rel_dist_d, lambda = post_rel_lambda_d, gamma = post_rel_gamma_d)
      g_mat[6, 8] <- matrix_vals$rate * hr_t_d
      s_mat[6, 8] <- exp(log(matrix_vals$surv) * hr_t_d)
    }
    
  } else if (dist_d != "llogistic") {  
    if (exists("hr_r_d_func")) {
      matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d * hr_t_d * hr_r_d_func(his_var), gamma = gamma_d)
    } else {
      matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d * hr_t_d * hr_r_d, gamma = gamma_d)
    }
    g_mat[6, 8] <- matrix_vals$rate
    s_mat[6, 8] <- matrix_vals$surv
    
  } else {
    matrix_vals <- assign_matrix_values(t = t, dist = dist_d, lambda = lambda_d, gamma = gamma_d)
    
    if (exists("hr_r_d_func")) {
      g_mat[6, 8] <- matrix_vals$rate * hr_t_d * hr_r_d_func(his_var)
      s_mat[6, 8] <- exp(log(matrix_vals$surv) * hr_t_d * hr_r_d_func(his_var))
      
    } else {
      g_mat[6, 8] <- matrix_vals$rate * hr_t_d * hr_r_d
      s_mat[6, 8] <- exp(log(matrix_vals$surv) * hr_t_d * hr_r_d) 
    }
  }
  
  diag(g_mat) <- -rowSums(g_mat)
  
  mat_list <- list(g_mat = g_mat, s_mat = s_mat)
  return(mat_list)
}

# Generate random holding time given that state_0 was arrived at time t_0
# Based on 'PINV' (Polynomial interpolation of INVerse CDF)
# Time scale: renewal
r_survival_renewal <- function(t_0, state_0, 
                               dists, lambdas, gammas, 
                               hr_t_r, hr_t_d, hr_r_d,
                               relapse_his = F, his_var = NA,
                               post_rel = F, post_rel_dist_d = NA, 
                               post_rel_lambda_d = NA, post_rel_gamma_d = NA) {
  # Simulation assumption
  if (state_0 %in% c(3, 4, 7, 8)) {
    stop("Absorbing state has been reached.")
  } 
  
  # Define cumulative distribution function given surv_t_0
  cdf <- function(t) {
    1 - prod(mat_list(t = t, dists = dists, lambdas = lambdas, gammas = gammas, 
                      hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d,
                      relapse_his = relapse_his, his_var = his_var,
                      post_rel = post_rel, post_rel_dist_d = post_rel_dist_d, 
                      post_rel_lambda_d = post_rel_lambda_d, post_rel_gamma_d = post_rel_gamma_d)[[2]][state_0, ])
  }
  
  # Set up random number generator
  gen <- pinv.new(cdf = cdf, lb = 0, ub = Inf, uresolution = 1e-12)
  
  # Generate random holding time
  holding_t <- uq(gen, runif(1, 0, 1))
  
  return(holding_t)
}

# Generate random holding time given that state_0 was arrived at time t_0
# Based on 'PINV' (Polynomial interpolation of INVerse CDF)
# Time scale: forward
r_survival_forward <- function(t_0, state_0, 
                               dists, lambdas, gammas, 
                               hr_t_r, hr_t_d, hr_r_d) {
  # Simulation assumption
  if (state_0 %in% c(3, 4, 7, 8)) {
    stop("Absorbing state has been reached.")
  } 
  
  if (all(dists == c("exponential", "exponential", "exponential"))) {
    rate_0 <- abs(mat_list(t = t_0, dists = dists, lambdas = lambdas, gammas = gammas, 
                        hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)[[1]][state_0, state_0])
    
    holding_t <- urexp(1, rate_0, lb = 0, ub = Inf)
  } else {
  # Conditional survival at t_0
  surv_t_0 <- mat_list(t = t_0, dists = dists, lambdas = lambdas, gammas = gammas, 
                       hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)[[2]][state_0, ]
  
  # Define cumulative distribution function given surv_t_0
  cdf <- function(t) {
    1 - prod(mat_list(t = t_0 + t, dists = dists, lambdas = lambdas, gammas = gammas, 
                      hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)[[2]][state_0, ] / surv_t_0)
  }
  
  # Set up random number generator
  gen <- pinv.new(cdf = cdf, lb = 0, ub = Inf, uresolution = 1e-12)
  
  # Generate random holding time
  holding_t <- uq(gen, runif(1, 0, 1))
  }
  return(holding_t)
}

# Return transition probability given the holding time (delta_t), the state arrived (state_0), and the time of this arrival (t_0)
# Time scale: renewal
trans_prob_renewal <- function(t_0, delta_t, state_0, 
                               dists, lambdas, gammas, 
                               hr_t_r, hr_t_d, hr_r_d,
                               relapse_his = F, his_var = NA,
                               post_rel = F, post_rel_dist_d = NA, 
                               post_rel_lambda_d = NA, post_rel_gamma_d = NA) {
  # Simulation assumption
  if (state_0 %in% c(3, 4, 7, 8)) {
    stop("Absorbing state has been reached.")
  }
  
  trans_prob <- mat_list(t = delta_t, dists = dists, lambdas = lambdas, gammas = gammas, 
                         hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d,
                         relapse_his = relapse_his, his_var = his_var,
                         post_rel = post_rel, post_rel_dist_d = post_rel_dist_d, 
                         post_rel_lambda_d = post_rel_lambda_d, post_rel_gamma_d = post_rel_gamma_d)[[1]][state_0, ]
  
  trans_prob[state_0] <- 0
  
  prob_sum <- sum(trans_prob)
  
  if (prob_sum > 0) {
    trans_prob <- trans_prob / prob_sum
  } else {
    stop("No valid transitions available. Sum of transition probabilities is zero.")
  }
  
  return(trans_prob)
}

# Return transition probability given the holding time (delta_t), the state arrived (state_0), and the time of this arrival (t_0)
# Time scale: forward
trans_prob_forward <- function(t_0, delta_t, state_0, 
                               dists, lambdas, gammas, 
                               hr_t_r, hr_t_d, hr_r_d) {
  # Simulation assumption
  if (state_0 %in% c(3, 4, 7, 8)) {
    stop("Absorbing state has been reached.")
  }
  
  trans_prob <- mat_list(t = t_0 + delta_t, dists = dists, lambdas = lambdas, gammas = gammas, 
                         hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)[[1]][state_0, ]
  
  trans_prob[state_0] <- 0
  
  prob_sum <- sum(trans_prob)
  
  if (prob_sum > 0) {
    trans_prob <- trans_prob / prob_sum
  } else {
    stop("No valid transitions available. Sum of transition probabilities is zero.")
  }
  
  return(trans_prob)
}

# Simulate one path
# state_0: initial status for the Markov jump process, set default is 1
# t_star: time up to which the Markov jump process is simulated
# Time scale: renewal
sim_path_renewal <- function(state_0 = 1, t_star, 
                             dists, lambdas, gammas, 
                             hr_t_r, hr_t_d, hr_r_d,
                             relapse_his = F, 
                             post_rel = F, post_rel_dist_d = NA, 
                             post_rel_lambda_d = NA, post_rel_gamma_d = NA) {
  # Initialize the path
  delta_t <- 0
  cum_t <- 0
  state_0 <- state_0
  res <- c(delta_t, cum_t, state_0)
  
  while (cum_t < t_star & state_0 %in% c(1, 2, 5, 6)) { # Ensure absorbing state has not been reached
    
    # Given state_0 was arrived at t_0 = cum_t
    t_0 <- cum_t
    
    # Assign value to his_var if relapse == T
    if (relapse_his) {
      his_var <- cum_t
    } else {his_var <- NA}
    
    # Step1: generate random holding time spent at state_0
    delta_t <- r_survival_renewal(t_0 = t_0, state_0 = state_0, 
                                  dists = dists, lambdas = lambdas, gammas = gammas, 
                                  hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d =  hr_r_d,
                                  relapse_his = relapse_his, his_var = his_var,
                                  post_rel = post_rel, post_rel_dist_d = post_rel_dist_d, 
                                  post_rel_lambda_d = post_rel_lambda_d, post_rel_gamma_d = post_rel_gamma_d)
    
    cum_t <- cum_t + delta_t
    
    # Step2: where to jump to
    trans_prob <- trans_prob_renewal(t_0 = t_0, delta_t = delta_t, state_0 = state_0, 
                                     dists = dists, lambdas = lambdas, gammas = gammas, 
                                     hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d =  hr_r_d,
                                     relapse_his = relapse_his, his_var = his_var,
                                     post_rel = post_rel, post_rel_dist_d = post_rel_dist_d, 
                                     post_rel_lambda_d = post_rel_lambda_d, post_rel_gamma_d = post_rel_gamma_d)
    
    if (any(trans_prob == 1)) {
      state_1 <- seq_along(trans_prob)[trans_prob == 1]
    } else {
      state_1 <- sample.int(length(trans_prob), 1, prob = trans_prob)
    }
    
    # Step3: trim the time and set state back if this jump is beyond t_star
    if (cum_t > t_star) {
      delta_t <- delta_t - (cum_t - t_star)
      cum_t <- t_star
      state_1 <- state_0
    }
    
    state_0 <- state_1
    res <- rbind(res, c(delta_t, cum_t, state_0))
  }
  
  res <- as.data.frame(res)
  colnames(res) <- c("delta_t", "total_t", "state")
  rownames(res) <- NULL
  
  return(res)
}

# Simulate one path
# state_0: initial status for the Markov jump process, set default is 1
# t_star: time up to which the Markov jump process is simulated
# Time scale: forward
sim_path_forward <- function(state_0 = 1, t_star, 
                             dists, lambdas, gammas, 
                             hr_t_r, hr_t_d, hr_r_d) {
  # Initialize the path
  delta_t <- 0
  cum_t <- 0
  state_0 <- state_0
  res <- c(delta_t, cum_t, state_0)
  
  while (cum_t < t_star & state_0 %in% c(1, 2, 5, 6)) { # Ensure absorbing state has not been reached
    
    # Given state_0 was arrived at t_0 = cum_t
    t_0 <- cum_t
    
    # Step1: generate random holding time spent at state_0
    delta_t <- r_survival_forward(t_0 = t_0, state_0 = state_0, 
                                  dists = dists, lambdas = lambdas, gammas = gammas, 
                                  hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d =  hr_r_d)
    cum_t <- cum_t + delta_t
    
    # Step2: where to jump to
    trans_prob <- trans_prob_forward(t_0 = t_0, delta_t = delta_t, state_0 = state_0, 
                                     dists = dists, lambdas = lambdas, gammas = gammas, 
                                     hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d =  hr_r_d)
    if (any(trans_prob == 1)) {
      state_1 <- seq_along(trans_prob)[trans_prob == 1]
    } else {
      state_1 <- sample.int(length(trans_prob), 1, prob = trans_prob)
    }
    
    # Step3: trim the time and set state back if this jump is beyond t_star
    if (cum_t > t_star) {
      delta_t <- delta_t - (cum_t - t_star)
      cum_t <- t_star
      state_1 <- state_0
    }
    
    state_0 <- state_1
    res <- rbind(res, c(delta_t, cum_t, state_0))
  }
  
  res <- as.data.frame(res)
  colnames(res) <- c("delta_t", "total_t", "state")
  rownames(res) <- NULL
  
  return(res)
}

# Prepare data set for time-dependent cox
create_cox_dat <- function(path_list, t_star) {
  
  # Create matrix for variables: treatment, relapse, and death
  var_mat <- cbind(T = c(rep(0, 4), rep(1, 4)), R = rep(c(0, 1), 4), D = rep(c(0, 0, 1, 1), 2))
  
  cox_dat <- NULL
  
  for (i in seq_along(path_list)) {
    path <- path_list[[i]]
    
    if (nrow(path) == 1) {
      df <- cbind(ID = i, Start = 0, Stop = t_star, Treat = 0, Relapse = 0, Death = 0)
      
    } else {
      Start <- path[1 : (nrow(path) - 1), 2]
      Stop <- path[2 : nrow(path), 2]
      TRD <- cbind(matrix(var_mat[path[-nrow(path), 3], 1 : 2], ncol = 2), var_mat[path[-1, 3], 3])
      df <- as.data.frame(cbind(i, Start, Stop, TRD))
      
      colnames(df) <- c("ID", "Start", "Stop", "Treat", "Relapse", "Death")
      rownames(df)<-NULL 
    }
    
    cox_dat <- rbind(cox_dat, df)
  }
  
  return(cox_dat)
}

# Prepare data set for supper landmark
# t_interval: length of the time interval
# t_horizon: time window for outcome to be predicted
create_landmark_dat <- function(path_list, t_star, t_interval, t_horizon) {
  
  # Initialize the data frame
  n_pat <- length(path_list)
  df <- as.data.frame(cbind(ID = seq_len(n_pat),
                            Treat_t = rep(t_star, n_pat), Treat_s = rep(0, n_pat),
                            Relapse_t = rep(t_star, n_pat), Relapse_s = rep(0, n_pat),
                            Death_t = rep(t_star, n_pat), Death_s = rep(0, n_pat)))
  
  for (i in seq_len(n_pat)) {
    path <- path_list[[i]]
    df_i <- df[df$ID == i, ]
    
    # first jump to "010"
    if (path[2, 3] == 2) {
      df_i$Relapse_t <- path[2, 2]
      df_i$Relapse_s <- 1
      
      if (nrow(path) > 2 & path[3, 3] == 4) {
        df_i$Death_t <- path[3, 2]
        df_i$Death_s <- 1
      }
      
      # first jump to "001"
    } else if (path[2, 3] == 3) {
      df_i$Death_t <- path[2, 2]
      df_i$Death_s <- 1
      
      # first jump to "100"
    } else if (path[2, 3] == 5) {
      df_i$Treat_t <- path[2, 2]
      df_i$Treat_s <- 1
      
      if (nrow(path) > 2 & path[3, 3] == 6) {
        df_i$Relapse_t <- path[3, 2]
        df_i$Relapse_s <- 1
        
        if (nrow(path) > 3 & path[4, 3] == 8) {
          df_i$Death_t <- path[4, 2]
          df_i$Death_s <- 1
        }
        
      } else if (nrow(path) > 2 & path[3, 3] == 7) {
        df_i$Death_t <- path[3, 2]
        df_i$Death_s <- 1
      }
    }
    df[df$ID == i, ] <- df_i
  }
  
  t_seq <- seq(0, t_star, by = t_interval)
  n_int <- length(t_seq)
  landmark_dat <- NULL
  
  # Sliding landmark for treatment
  for (i in 1 : (n_int - 1)) {
    landmark_dat <- rbind(landmark_dat, cutLM(data = df, outcome = list(time = "Death_t", status = "Death_s"),
                                              LM = t_seq[i], horizon = min(t_seq[i] + t_horizon, t_star),
                                              covs = list(fixed = "ID", varying = "Treat_t")))
  }
  colnames(landmark_dat)[which(colnames(landmark_dat) == "Treat_t")] <- "Treat"
  
  # Add relapse
  landmark_dat$Relapse <- rep(NA, nrow(landmark_dat))
  for (i in seq_len(n_pat)) {
    Relapse_i <- df$Relapse_t[df$ID == i]
    landmark_dat$Relapse[landmark_dat$ID == i] <- as.numeric(landmark_dat$LM[landmark_dat$ID == i] >= Relapse_i)
  }
  
  landmark_dat$Landmark <- landmark_dat$LM
  return(landmark_dat)
}

# Prepare data set for gformula
# t_interval: length of the time interval
create_gform_dat <- function(path_list, t_star, t_interval) {
  
  # Initialize the data frame
  t_seq <- seq(0, t_star, by = t_interval)
  n_int <- length(t_seq)
  n_pat <- length(path_list)
  
  df <- as.data.frame(cbind(ID = rep(seq_len(n_pat), each = n_int),
                            Time = rep(t_seq, n_pat),
                            Treat = rep(0, n_int * n_pat),
                            Relapse = rep(0, n_int * n_pat),
                            Death = rep(0, n_int * n_pat)))
  
  for (i in seq_len(n_pat)) {
    
    path <- path_list[[i]]
    df_i <- df[df$ID == i, ]
    
    # first jump to "010"
    if (path[2, 3] == 2) {
      df_i$Relapse[(ceiling(path[2, 2] / t_interval) + 1) : n_int] <- 1
      if (nrow(path) > 2 & path[3, 3] == 4) {
        df_i$Death[(ceiling(path[3, 2] / t_interval) + 1) : n_int] <- 1
      }
      
      # first jump to "001"
    } else if (path[2, 3] == 3) {
      df_i$Death[(ceiling(path[2, 2] / t_interval) + 1) : n_int] <- 1
      
      # first jump to "100"
    } else if (path[2, 3] == 5) {
      df_i$Treat[(ceiling(path[2, 2] / t_interval) + 1) : n_int] <- 1
      if (nrow(path) > 2 & path[3, 3] == 6) {
        df_i$Relapse[(ceiling(path[3, 2] / t_interval) + 1) : n_int] <- 1
        if (nrow(path) > 3 & path[4, 3] == 8) {
          df_i$Death[(ceiling(path[4, 2] / t_interval) + 1) : n_int] <- 1
        } 
      } else if (nrow(path) > 2 & path[3, 3] == 7) {
        df_i$Death[(ceiling(path[3, 2] / t_interval) + 1) : n_int] <- 1
      }
    }
    df[df$ID == i, ] <- df_i
  }
  
  gform_dat <- df %>% group_by(ID) %>% mutate(ind = cumsum(Death)) %>% 
    filter(ind %in% 0:1) %>% select(-ind)
  
  return(gform_dat)
}

# Hypothetical history for time-dependent cox prediction
create_cox_int <- function(sim_dat_int, int_time = int_time, t_star = t_star) {
  int_row <- data.frame(Start = 0, Stop  = 0, Treat = 0, Relapse = 0, Var1 = 0, Var2 = 0)
  dat_cox_int <- NULL
  
  if (int_time >= t_star) {
    for (i in seq_len(nrow(sim_dat_int))) {
      if (sim_dat_int$Relapse_t[i] < t_star) {
        dat_cox_int_i <- int_row[c(1, 1), ]
        dat_cox_int_i$Relapse[2] <- 1
        dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- sim_dat_int$Relapse_t[i]
        dat_cox_int_i$Stop[2] <- t_star
      } else {
        dat_cox_int_i <- int_row
        dat_cox_int_i$Stop <- t_star
      }
      dat_cox_int_i$Var1 <- sim_dat_int$Var1[i]
      dat_cox_int_i$Var2 <- sim_dat_int$Var2[i]
      dat_cox_int_i$ID <- i
      
      dat_cox_int <- rbind(dat_cox_int, dat_cox_int_i)
    }
  } else {
    for (i in seq_len(nrow(sim_dat_int))) {
      if (sim_dat_int$Treat_s[i] == 1) {
        if (sim_dat_int$Relapse_t[i] < int_time) {
          dat_cox_int_i <- int_row[c(1, 1, 1), ]
          dat_cox_int_i$Relapse[2:3] <- 1
          dat_cox_int_i$Treat[3] <- 1
          dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- sim_dat_int$Relapse_t[i]
          dat_cox_int_i$Start[3] <- dat_cox_int_i$Stop[2] <- int_time
          dat_cox_int_i$Stop[3] <- t_star
        } else if (sim_dat_int$Relapse_t[i] == int_time) {
          dat_cox_int_i <- int_row[c(1, 1), ]
          dat_cox_int_i$Relapse[2] <- 1
          dat_cox_int_i$Treat[2] <- 1
          dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- sim_dat_int$Relapse_t[i]
          dat_cox_int_i$Stop[2] <- t_star
        } else if (sim_dat_int$Relapse_t[i] < t_star) {
          dat_cox_int_i <- int_row[c(1, 1, 1), ]
          dat_cox_int_i$Relapse[3] <- 1
          dat_cox_int_i$Treat[2:3] <- 1
          dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- int_time
          dat_cox_int_i$Start[3] <- dat_cox_int_i$Stop[2] <- sim_dat_int$Relapse_t[i]
          dat_cox_int_i$Stop[3] <- t_star
        } else {
          dat_cox_int_i <- int_row[c(1, 1), ]
          dat_cox_int_i$Treat[2] <- 1
          dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- int_time
          dat_cox_int_i$Stop[2] <- t_star
        }
        dat_cox_int_i$Var1 <- sim_dat_int$Var1[i]
        dat_cox_int_i$Var2 <- sim_dat_int$Var2[i]
        dat_cox_int_i$ID <- i
        
        dat_cox_int <- rbind(dat_cox_int, dat_cox_int_i)
      } else {
        if (sim_dat_int$Relapse_t[i] < t_star) {
          dat_cox_int_i <- int_row[c(1, 1), ]
          dat_cox_int_i$Relapse[2] <- 1
          dat_cox_int_i$Start[2] <- dat_cox_int_i$Stop[1] <- sim_dat_int$Relapse_t[i]
          dat_cox_int_i$Stop[2] <- t_star
        } else {
          dat_cox_int_i <- int_row
          dat_cox_int_i$Stop <- t_star
        }
        dat_cox_int_i$Var1 <- sim_dat_int$Var1[i]
        dat_cox_int_i$Var2 <- sim_dat_int$Var2[i]
        dat_cox_int_i$ID <- i
        
        dat_cox_int <- rbind(dat_cox_int, dat_cox_int_i)
      }
    }
  }
  return(dat_cox_int)
}

# Hypothetical history for landmark prediction
create_landmark_int <- function(sim_dat_int, int_time = int_time, t_star = t_star) {
  
  int_row <- data.frame(Treat = rep(0, t_star), Relapse = rep(0, t_star), 
                        LM = 0 : 35, Landmark = 0 : 35, 
                        Var1 = rep(0, t_star), Var2 = rep(0, t_star))
  
  dat_landmark_int <- NULL
  
  if (int_time > (t_star - 1)) {
    for (i in seq_len(nrow(sim_dat_int))) {
      dat_landmark_int_i <- int_row
      
      if (sim_dat_int$Relapse_t[i] <= (t_star - 1)) {
        dat_landmark_int_i$Relapse[(ceiling(sim_dat_int$Relapse_t[i]) + 1) : t_star] <- 1
      } 
      
      dat_landmark_int_i$Var1 <- sim_dat_int$Var1[i]
      dat_landmark_int_i$Var2 <- sim_dat_int$Var2[i]
      dat_landmark_int_i$ID <- i
      
      dat_landmark_int <- rbind(dat_landmark_int, dat_landmark_int_i)
    }
  } else {
    for (i in seq_len(nrow(sim_dat_int))) {
      dat_landmark_int_i <- int_row
      
      if (sim_dat_int$Treat_s[i] == 1) {
        dat_landmark_int_i$Treat[(ceiling(int_time) + 1) : t_star] <- 1
      }
      
      if (sim_dat_int$Relapse_t[i] <= (t_star - 1)) {
        dat_landmark_int_i$Relapse[(ceiling(sim_dat_int$Relapse_t[i]) + 1) : t_star] <- 1
      } 
      
      dat_landmark_int_i$Var1 <- sim_dat_int$Var1[i]
      dat_landmark_int_i$Var2 <- sim_dat_int$Var2[i]
      dat_landmark_int_i$ID <- i
      
      dat_landmark_int <- rbind(dat_landmark_int, dat_landmark_int_i)
    }
  }
  return(dat_landmark_int)
}

# Prepare wide format data
create_wide_dat <- function(path_list, t_star) {
  
  # Initialize the data frame
  n_pat <- length(path_list)
  df <- as.data.frame(cbind(ID = seq_len(n_pat),
                            Treat_t = rep(t_star, n_pat), Treat_s = rep(0, n_pat),
                            Relapse_t = rep(t_star, n_pat), Relapse_s = rep(0, n_pat),
                            Death_t = rep(t_star, n_pat), Death_s = rep(0, n_pat)))
  
  for (i in seq_len(n_pat)) {
    path <- path_list[[i]]
    df_i <- df[df$ID == i, ]
    
    # first jump to "010"
    if (path[2, 3] == 2) {
      df_i$Relapse_t <- path[2, 2]
      df_i$Relapse_s <- 1
      if (nrow(path) > 2 & path[3, 3] == 4) {
        df_i$Death_t <- path[3, 2]
        df_i$Death_s <- 1
      }
      
      # first jump to "001"
    } else if (path[2, 3] == 3) {
      df_i$Death_t <- path[2, 2]
      df_i$Death_s <- 1
      
      # first jump to "100"
    } else if (path[2, 3] == 5) {
      df_i$Treat_t <- path[2, 2]
      df_i$Treat_s <- 1
      if (nrow(path) > 2 & path[3, 3] == 6) {
        df_i$Relapse_t <- path[3, 2]
        df_i$Relapse_s <- 1
        if (nrow(path) > 3 & path[4, 3] == 8) {
          df_i$Death_t <- path[4, 2]
          df_i$Death_s <- 1
        } 
      } else if (nrow(path) > 2 & path[3, 3] == 7) {
        df_i$Death_t <- path[3, 2]
        df_i$Death_s <- 1
      }
    }
    df_i$Treat_t <- min(df_i$Treat_t, df_i$Relapse_t, df_i$Death_t)
    df_i$Relapse_t <- min(df_i$Relapse_t, df_i$Death_t)
    df[df$ID == i, ] <- df_i
  }
  
  return(df)
}

# Prepare wide format data for multi-state modelling to obtain Aalen-Johansen estimator
create_wide_mstate_dat <- function(path_list, t_star) {
  
  # Initialize the data frame
  n_pat <- length(path_list)
  df <- as.data.frame(cbind(ID = seq_len(n_pat),
                            t_2 = rep(t_star, n_pat), s_2 = rep(0, n_pat),
                            t_3 = rep(t_star, n_pat), s_3 = rep(0, n_pat),
                            t_4 = rep(t_star, n_pat), s_4 = rep(0, n_pat),
                            t_5 = rep(t_star, n_pat), s_5 = rep(0, n_pat),
                            t_6 = rep(t_star, n_pat), s_6 = rep(0, n_pat),
                            t_7 = rep(t_star, n_pat), s_7 = rep(0, n_pat),
                            t_8 = rep(t_star, n_pat), s_8 = rep(0, n_pat)))
  
  for (i in 1:n_pat) {
    path <- path_list[[i]]
    df_i <- df[df$ID == i, ]
    
    if (path$state[length(path$state)] == path$state[length(path$state) - 1]) {
      path <- path[-c(1, length(path$state)), ]
    } else {
      path <- path[-1, ]
    }
    
    for (s in path$state) {
      df_i[paste0("s_", s)] <- 1
      df_i[paste0("t_", s)] <- path$total_t[path$state == s]
    }
    df[df$ID == i, ] <- df_i
  }
  return(df)
}





