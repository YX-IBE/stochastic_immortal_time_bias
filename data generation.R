# S1: Weibull distribution, no treatment effect
# Time scale: forward

setwd("~/S1")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                        hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                        hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                        hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                    hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# S2: Weibull distribution, treatment effect only on death
# Time scale: forward

setwd("~/S2")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# S3: Weibull distribution, treatment effect on both relapse and death
# Time scale: forward

setwd("~/S3")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "weibull")
lambdas_1 <- c(0.012, 0.12, 0.048)
lambdas_2 <- c(0.012, 0.16, 0.060)
lambdas_3 <- c(0.012, 0.12, 0.072)
lambdas_4 <- c(0.012, 0.16, 0.084)
gammas <- c(1.5, 0.55, 0.60)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# S4: Log-logistic distribution, no treatment effect
# Time scale: forward

setwd("~/S4")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 1
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# S5: Log-logistic distribution, treatment effect only on death
# Time scale: forward

setwd("~/S5")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 1
hr_t_d <- 0.5
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# S6: Log-logistic distribution, treatment effect on both relapse and death
# Time scale: forward

setwd("~/S6")
Sys.setenv(OPENBLAS_NUM_THREADS = 1)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

library(doParallel)
library(foreach)
library(doRNG)

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/forward")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_forward(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/forward/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()

# Time scale: renewal

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- 2

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/renewal")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/renewal/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()


# Time scale: extended

# Parameters for Markov process
n_pat <- 250
t_star <- 36
dists <- c("weibull", "weibull", "llogistic")
lambdas_1 <- c(0.02, 0.12, 0.012)
lambdas_2 <- c(0.02, 0.16, 0.0125)
lambdas_3 <- c(0.02, 0.12, 0.013)
lambdas_4 <- c(0.02, 0.16, 0.0135)
gammas <- c(1.5, 0.55, 1.3)
hr_t_r <- 0.5
hr_t_d <- 0.5
hr_r_d <- "2 * exp(log(0.5 + 0.5 / (1 + exp((his_var - 6)))))"
relapse_his <- T

# Parameters for parallel simulations
n_simul <- 1000
n_batch <- 10
batch_size <- n_simul / n_batch

n_cores <- 112
cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterEvalQ(cl, source("./code/functions.R"))

RNG <- 1809

dirs <- c("./simulation", 
          "./simulation/extended")

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir)
    print(paste("Directory", dir, "created."))
  } else {
    print(paste("Directory", dir, "already exists."))
  }
}

for (batch in 1:n_batch) {
  print(paste("Processing batch", batch))
  batch_rng <- RNG + batch * 100
  
  batch_time <- system.time({
    path_lists_batch <- foreach(simul = ((batch - 1) * batch_size + 1):(batch * batch_size),
                                .combine = 'list', .multicombine = TRUE, .options.RNG = batch_rng) %dorng% {
                                  
                                  path_list <- vector("list", n_pat * 4)
                                  
                                  for (id in 1:n_pat) {
                                    path_list[[id * 4 - 3]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_1, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 2]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_2, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4 - 1]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_3, gammas = gammas,
                                                                                hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                    path_list[[id * 4]] <- sim_path_renewal(t_star = t_star, dists = dists, lambdas = lambdas_4, gammas = gammas,
                                                                            hr_t_r = hr_t_r, hr_t_d = hr_t_d, hr_r_d = hr_r_d, relapse_his = relapse_his)
                                  }
                                  return(path_list)
                                }
    saveRDS(path_lists_batch, file = paste0("./simulation/extended/path_lists_batch_", batch, ".rds"), compress = TRUE)
  })
  
  print(paste("Finished batch", batch, "in", batch_time[3] / 60, "minutes"))
}

stopCluster(cl)
registerDoSEQ()



