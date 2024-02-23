library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

Table_tau <- function(n, p, alpha) {
  Rep <- 10000
  
  PMLE <- matrix(ncol = Rep, nrow = length(p) - 1)
  
  for (i in 1:Rep) {
    fit_tau <- 1
    while (all(fit_tau == 1)) {
      x <- rpwexp(n, p, alpha)
      
      x_min <- min(x)
      x_max <- max(x)
      d1 <- 3.6
      d2 <- 4.4
      A_2 <- matrix(c(1, -1, 1, 0, 0, 1, -1, -1), ncol = 2)
      d_2 <- c(-x_min, -d1, d2, x_max)
      
      profile <- function(y) profile_loglik_pwexp_bc(p = c(p[1], y), x)
      
      fit_tau <- tryCatch(
        maxSANN(
          profile,
          start = c(x_min + 0.05, x_min + d1 + 0.1),
          constraints = list(ineqA = A_2, ineqB = d_2),
          control = list(iterlim = 500)
        )$estimate, error = function (t) return(1)
      )
      
      fit_tau
    }
    
    PMLE[, i] <- fit_tau
  }
  
  BIAS <- apply(PMLE - p[-1], 1, mean)
  MSE <- apply(PMLE - p[-1], 1, function(x) mean(x^2))
  
  data.frame("n" = paste(n), BIAS, MSE)
}

# Scenario 3 --------------------------------------------------------------

set.seed(2024)
number_of_cores <- 8
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 8)

ptm <- proc.time() #start

estimates <- 
  foreach(n_grid = seq(10, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_tau(n_grid, p = c(4, 6, 10), alpha = c(0.3, 0.5, 1.9)) 

proc.time() - ptm #final

temp <- purrr::map_dfr(estimates, rbind)[
  c(seq(1, 40, by = 2), seq(1, 40, by = 2) + 1),]

print.xtable(xtable(temp, digits = 3), include.rownames = F)
