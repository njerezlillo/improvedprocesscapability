library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

set.seed(2024)

Table_CpmA <- function(n, p, alpha, USL, LSL, Tc) {
  Rep <- 10000
  
  M <- qpwexp(0.5, p, alpha)
  Lp <- qpwexp(0.00135, p, alpha)
  Up <- qpwexp(0.99865, p, alpha)
  temp <- min(USL - Tc, Tc - LSL)
  CpmA_real <- temp/(3 * sqrt(((Up - Lp)/3)^2 + (M - Tc)^2)) 
  
  Upc <- Lpc <- Mc <- CpmA <- EE <- matrix(ncol = Rep, nrow = 1)
  
  for (i in 1:Rep) {
    x <- rpwexp(n, p, alpha)
    est <- mle_pwexp_bc(x, p)
    H <- diag(est^2 / (n * diff_c(auxiliar_pwexp(p, est))))
    Upc[, i] <- qpwexp(0.99865, p, est)
    Lpc[, i] <- qpwexp(0.00135, p, est)
    Mc[, i] <- qpwexp(0.5, p, est)
    aux <- min(USL - Tc, Tc - LSL)
    CpmA[, i] <- 
      aux/(3 * sqrt(((Upc[, i] - Lpc[, i])/3)^2 + (Mc[, i] - Tc)^2))
    EE[, i] <- md_CpmA(est, H, p, LSL, USL, Tc)
  }
  
  li <- CpmA - 1.96 * sqrt(EE)
  ls <- CpmA + 1.96 * sqrt(EE)
  
  MRE <- mean(CpmA / CpmA_real)
  MSE <- mean((CpmA - CpmA_real)^2)
  CP <- apply(li < CpmA_real & ls > CpmA_real, 1, mean)
  
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
  Table_CpmA(n_grid, p = c(2, 5, 10), alpha = c(0.1, 0.05, 10), USL = 15, LSL = 0.1, Tc = 7.93) 

proc.time() - ptm #final

xtable(data.frame(t(sapply(estimates, c))), digits = 3)