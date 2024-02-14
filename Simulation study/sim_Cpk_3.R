library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

set.seed(2024)

Table_Cpk <- function(n, p, alpha, USL, LSL) {
  Rep <- 10000
  
  M <- qpwexp(0.5, p, alpha)
  Lp <- qpwexp(0.00135, p, alpha)
  Up <- qpwexp(0.99865, p, alpha)
  Cpk_real <- min((USL - M)/(Up - M), (M - LSL)/(M - Lp))
  Cpk_real
  
  Upc <- Lpc <- Mc <- Cpk <- EE <- matrix(ncol = Rep, nrow = 1)
  
  i <- 1
  while (i <= Rep) {
    x <- rpwexp(n, p, alpha)
    est <- mle_pwexp_bc(x, p)
    if (all(est > 0)) {
      H <- diag(est^2 / (n * diff_c(auxiliar_pwexp(p, est))))
      Upc[, i] <- qpwexp(0.99865, p, est)
      Lpc[, i] <- qpwexp(0.00135, p, est)
      Mc[, i] <- qpwexp(0.5, p, est)
      Cpk[, i] <- min((USL - Mc[, i])/(Upc[, i] - Mc[, i]),
                    (Mc[, i] - LSL)/(Mc[, i] - Lpc[, i]))
      EE[, i] <- md_Cpk(est, H, p, LSL, USL)
      i <- i+1
    }
  }
  
  li <- Cpk - 1.96 * sqrt(EE)
  ls <- Cpk + 1.96 * sqrt(EE)
  
  MRE <- mean(Cpk / Cpk_real)
  MSE <- mean((Cpk - Cpk_real)^2)
  CP <- apply(li < Cpk_real & ls > Cpk_real, 1, mean)
  
  data.frame("n" = paste(n), "MRE" = MRE, "MSE" = MSE, "CP" = CP)
}

number_of_cores <- 20 # 40
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 20) # 40

ptm <- proc.time() #start

estimates <- 
  foreach(n_grid = seq(20, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_Cpk(n_grid, p = c(4, 6, 10), alpha = c(0.3, 0.5, 1.9), USL = 15, LSL = 0.1) 

proc.time() - ptm #final

xtable(data.frame(t(sapply(estimates, c))), digits = 3)