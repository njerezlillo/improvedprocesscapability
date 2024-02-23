library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

set.seed(2024)

Table_alpha <- function(n, p, alpha) {
  Rep <- 10000
  
  MLE <- matrix(ncol = Rep, nrow = (length(p)))
  MLE_bc <- matrix(ncol = Rep, nrow = (length(p)))
  
  for (i in 1:Rep) {
    while (TRUE) {
      x <- rpwexp(n, p, alpha)
      est <- mle_pwexp_bc(x, p)
      if(all(est > 0)) break 
    }
    
    MLE[, i] <- mle_pwexp(x, p)
    MLE_bc[, i] <- mle_pwexp_bc(x, p)  
  }
  
  BIAS <- apply((MLE - alpha), 1, mean)
  MSE <- apply((MLE - alpha), 1, function(x) mean(x^2))
  
  BIAS_bc <- apply((MLE_bc - alpha), 1, mean)
  MSE_bc <- apply((MLE_bc - alpha), 1, function(x) mean(x^2))
  
  data.frame("n" = paste(n), 
             BIAS = paste0(sprintf(BIAS_bc, fmt = '%.3f'), " (", sprintf(BIAS, fmt = '%.3f'), ")"), 
             MSE = paste0(sprintf(MSE_bc, fmt = '%.3f'), " (", sprintf(MSE, fmt = '%.3f'), ")"))
}

# Scenario 1 --------------------------------------------------------------

set.seed(2024)
number_of_cores <- 8
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 8)

ptm <- proc.time() #start

estimates_1 <- 
  foreach(n_grid = seq(10, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_alpha(n_grid, p = c(1, 2, 3), alpha = c(0.4, 0.6, 0.8)) 

proc.time() - ptm #final

# Scenario 2 --------------------------------------------------------------

set.seed(2024)
number_of_cores <- 8
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 8)

ptm <- proc.time() #start

estimates_2 <- 
  foreach(n_grid = seq(10, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_alpha(n_grid, p = c(2, 5, 10), alpha = c(0.1, 0.05, 10)) 

proc.time() - ptm #final

# Scenario 3 --------------------------------------------------------------

set.seed(2024)
number_of_cores <- 8
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 8)

ptm <- proc.time() #start

estimates_3 <- 
  foreach(n_grid = seq(10, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_alpha(n_grid, p = c(4, 6, 10), alpha = c(0.3, 0.5, 1.9))

proc.time() - ptm #final

# Out ---------------------------------------------------------------------

temp <- cbind(
  data.frame(rep(NA, 20), rep(NA, 20)),
  purrr::map_dfr(estimates_1, rbind),
  data.frame(rep(NA, 20)),
  purrr::map_dfr(estimates_2, rbind)[,-1],
  data.frame(rep(NA, 20)),
  purrr::map_dfr(estimates_3, rbind)[,-1]
)

print.xtable(xtable(temp[c(seq(1, 60, by = 3),
                           seq(1, 60, by = 3) + 1,
                           seq(1, 60, by = 3) + 2), ]),
             include.rownames = F)
