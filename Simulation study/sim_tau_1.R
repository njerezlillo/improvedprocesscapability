library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("pwexp.R")
source("capability.R")

set.seed(2024)

Table_tau <- function(n, p, alpha) {
  Rep <- 10000
  
  PMLE <- matrix(ncol = Rep, nrow = length(p) - 1)
  
  for (i in 1:Rep) {
    fit_tau <- 1
    while (all(fit_tau == 1)) {
      x <- rpwexp(n, p, alpha)
      
      x_min <- quantile(sort(x), 0.3) # min(x)
      x_max <- quantile(sort(x), 0.7)
      d <- 0.8
      A_2 <- matrix(c(1, -1, 0, 0, 1, -1), ncol = 2)
      d_2 <- c(-x_min, -d, x_max)
      
      profile <- function(y) profile_loglik_pwexp(p = c(p[1], y), x)
      
      
      fit_tau <- tryCatch(
        maxSANN(
          profile,
          start = sort(runif(2, x_min, x_max)),
          constraints = list(ineqA = A_2, ineqB = d_2)
        )$estimate, error = function (t) return(1)
      )
    }
    
    PMLE[, i] <- fit_tau
  }
  
  BIAS <- apply(PMLE - p[-1], 1, mean)
  MSE <- apply(PMLE - p[-1], 1, function(x) mean(x^2))
  
  data.frame("n" = paste(n), BIAS, MSE)
}

number_of_cores <- 40 # 40
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 40) # 40

ptm <- proc.time() #start

estimates <- 
  foreach(n_grid = seq(20, 200, by = 10), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_tau(n_grid, alpha = c(0.4, 0.6, 0.8), p = c(1, 2, 3)) 

proc.time() - ptm #final

# beepr::beep(8)

out_1 <- rbind(
  data.frame(estimates[[1]][1, ]),
  data.frame(estimates[[2]][1, ]),
  data.frame(estimates[[3]][1, ]),
  data.frame(estimates[[4]][1, ]),
  data.frame(estimates[[5]][1, ]),
  data.frame(estimates[[6]][1, ]),
  data.frame(estimates[[7]][1, ]),
  data.frame(estimates[[8]][1, ]),
  data.frame(estimates[[9]][1, ]),
  data.frame(estimates[[10]][1, ]),
  data.frame(estimates[[11]][1, ]),
  data.frame(estimates[[12]][1, ]),
  data.frame(estimates[[13]][1, ]),
  data.frame(estimates[[14]][1, ]),
  data.frame(estimates[[15]][1, ]),
  data.frame(estimates[[16]][1, ]),
  data.frame(estimates[[17]][1, ]),
  data.frame(estimates[[18]][1, ]),
  data.frame(estimates[[19]][1, ])
) 

out_2 <- rbind(
  data.frame(estimates[[1]][2, ]),
  data.frame(estimates[[2]][2, ]),
  data.frame(estimates[[3]][2, ]),
  data.frame(estimates[[4]][2, ]),
  data.frame(estimates[[5]][2, ]),
  data.frame(estimates[[6]][2, ]),
  data.frame(estimates[[7]][2, ]),
  data.frame(estimates[[8]][2, ]),
  data.frame(estimates[[9]][2, ]),
  data.frame(estimates[[10]][2, ]),
  data.frame(estimates[[11]][2, ]),
  data.frame(estimates[[12]][2, ]),
  data.frame(estimates[[13]][2, ]),
  data.frame(estimates[[14]][2, ]),
  data.frame(estimates[[15]][2, ]),
  data.frame(estimates[[16]][2, ]),
  data.frame(estimates[[17]][2, ]),
  data.frame(estimates[[18]][2, ]),
  data.frame(estimates[[19]][2, ])
)

xtable(rbind(out_1, out_2), digits = 3)
