library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwexp.R")
source("capability.R")

Table_Cpk <- function(n, p, alpha, USL, LSL) {
  Rep <- 10000
  
  M <- qpwexp(0.5, p, alpha)
  Lp <- qpwexp(0.00135, p, alpha)
  Up <- qpwexp(0.99865, p, alpha)
  Cpk_real <- min((USL - M)/(Up - M), (M - LSL)/(M - Lp))
  Cpk_real
  
  Upc <- Lpc <- Mc <- Cpk <- EE <- matrix(ncol = Rep, nrow = 1)
  
  for (i in 1:Rep) {
    while (TRUE) {
      x <- rpwexp(n, p, alpha)
      est <- mle_pwexp_bc(x, p)
      if(all(est > 0)) break 
    }
    H <- diag(est^2 / (n * diff_c(auxiliar_pwexp(p, est))))
    Upc[, i] <- qpwexp(0.99865, p, est)
    Lpc[, i] <- qpwexp(0.00135, p, est)
    Mc[, i] <- qpwexp(0.5, p, est)
    Cpk[, i] <- min((USL - Mc[, i])/(Upc[, i] - Mc[, i]),
                    (Mc[, i] - LSL)/(Mc[, i] - Lpc[, i]))
    EE[, i] <- md_Cpk(est, H, p, LSL, USL)
  }
  
  li <- Cpk - 1.96 * sqrt(EE)
  ls <- Cpk + 1.96 * sqrt(EE)
  
  BIAS <- mean((Cpk - Cpk_real))
  MSE <- mean((Cpk - Cpk_real)^2)
  CP <- apply(li < Cpk_real & ls > Cpk_real, 1, mean)
  
  data.frame("n" = paste(n), BIAS, MSE, CP)
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
  Table_Cpk(n_grid, p = c(1, 2, 3), alpha = c(0.4, 0.6, 0.8), USL = 7, LSL = 1.3) 

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
  Table_Cpk(n_grid, p = c(2, 5, 10), alpha = c(0.1, 0.05, 10), USL = 15, LSL = 0.1) 

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
  Table_Cpk(n_grid, p = c(4, 6, 10), alpha = c(0.3, 0.5, 1.9), USL = 15, LSL = 0.1)

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

print.xtable(xtable(temp, digits = 3), include.rownames = F)
