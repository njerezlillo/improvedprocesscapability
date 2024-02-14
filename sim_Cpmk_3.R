library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

set.seed(2024)

Table_Cpmk <- function(n, p, alpha, USL, LSL, Tc) {
  Rep <- 10000
  
  M <- qpwexp(0.5, p, alpha)
  Lp <- qpwexp(0.00135, p, alpha)
  Up <- qpwexp(0.99865, p, alpha)
  temp1 <- (USL - M)/(3 * sqrt(((Up - M)/3)^2 + (M - Tc)^2)) 
  temp2 <- (M - LSL)/(3 * sqrt(((M - Lp)/3)^2 + (M - Tc)^2))
  Cpmk_real <- min(temp1, temp2)
  
  Upc <- Lpc <- Mc <- Cpmk <- EE <- matrix(ncol = Rep, nrow = 1)
  
  i <- 1
  while (i <= Rep) {
    x <- rpwexp(n, p, alpha)
    est <- mle_pwexp_bc(x, p)
    if (all(est > 0)) {
      H <- diag(est^2 / (n * diff_c(auxiliar_pwexp(p, est))))
      Upc[, i] <- qpwexp(0.99865, p, est)
      Lpc[, i] <- qpwexp(0.00135, p, est)
      Mc[, i] <- qpwexp(0.5, p, est)
      aux1 <- 
        (USL - Mc[, i])/(3 * sqrt(((Upc[, i] - Mc[, i])/3)^2 + (Mc[, i] - Tc)^2))
      aux2 <- 
        (Mc[, i] - LSL)/(3 * sqrt(((Mc[, i] - Lpc[, i])/3)^2 + (Mc[, i] - Tc)^2))
      Cpmk[, i] <- min(aux1, aux2)
      EE[, i] <- md_Cpmk(est, H, p, LSL, USL, Tc)
      i <- i + 1
    }      
  }
  
  li <- Cpmk - 1.96 * sqrt(EE)
  ls <- Cpmk + 1.96 * sqrt(EE)
  
  MRE <- mean(Cpmk / Cpmk_real)
  MSE <- mean((Cpmk - Cpmk_real)^2)
  CP <- apply(li < Cpmk_real & ls > Cpmk_real, 1, mean)
  
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
  Table_Cpmk(n_grid, p = c(4, 6, 10), alpha = c(0.1, 0.5, 1.9), USL = 15, LSL = 0.1, Tc = 7.28) 

proc.time() - ptm #final

xtable(data.frame(t(sapply(estimates, c))), digits = 3)