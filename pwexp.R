# C's
auxiliar_pwexp <- function (p, alpha)
{
  if (any(sort(p) != p)) {
    stop(paste("p must be crecient", "\n", ""))
  }
  if (any(alpha < 0)) {
    stop(paste("alpha must be > 0", "\n", ""))
  }
  
  aux <- c(1, rep(NA, length(p) - 1), 0)
  for (i in 2:length(p)) {
    aux[i] <- prod(exp(-alpha[i - 1]*(p[i] - p[i - 1])), aux[(i - 1)])
  }
  
  return (aux)
}

# density distribution function
dpwexp <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwexp(p, alpha)
  j <- max(which(p < x))
  d <- alpha[j] * exp(-alpha[j] * (x - p[j])) * C[j]
  
  return(d)
}

# cumulative distribution function
ppwexp <- function (x, p, alpha)
{
  C <- auxiliar_pwexp(p, alpha)
  if (x < p[1]) { 
    u = 0 
  } else {
    j <- max(which(p <= x))
    u <- 1 - exp(-alpha[j] * (x - p[j])) * C[j]
  }
  return(u)
}

# survival function
spwexp <- function (x, p, alpha)
{
  s <- 1 - ppwexp(x, p, alpha)
  return(s)
}

# hazard function
hpwexp <- function (x, p, alpha)
{
  if(x<=p[1]){ stop(return(0)) }
  C <- auxiliar_pwexp(p, alpha)
  j <- max(which(p <= x))
  alpha[j]
}

# quantile function
qpwexp <- function (u, p, alpha)
{
  C <- auxiliar_pwexp(p, alpha)
  aux <- u < (1 - C)
  j <- min(which(aux == TRUE)) - 1
  x <- p[j] - log((1 - u)/ C[j]) / alpha[j]
  
  return(x)
}

# identifying the observations in each interval
index_each_interval <- function (x, p)
{
  k <- length(p)
  nj <- vector("list", length = k)
  for (j in 1:k)
  {
    if (j < k)
    {
      nj[[j]] <- which(x >= p[j] & x < p[j + 1])
    } else {
      nj[[j]] <- which(x >= p[j])
    }
  }
  
  return (nj)
}

# counting the observations in each interval
n_each_interval <- function (x, p)
{
  nj <- lengths(index_each_interval(x, p))
  return (nj)
}

# sampling piece-wise exponential
rpwexp <- function(n, p, alpha)
{
  q <- Vectorize(function(x) qpwexp(x, p, alpha), "x")
  while(TRUE)
  {
    x <- q(runif(n))
    nj <- n_each_interval(x, p)
    if (all(nj > 2)) break
  }
  
  return (x)
}

# Simulation under random censoring
censrand_pwexp <- function(n, p, alpha, PI)
{
  f <- function(x) x * dpwexp(x, p, alpha)
  E <- integrate(Vectorize(f), p[1], Inf)$value
  lambda <- E/PI
  q <- Vectorize(function(x) qpwexp(x, p, alpha), "x")
  
  while(TRUE)
  {
    x <- q(runif(n))
    v <- runif(n, 0, lambda)
    delta <- (x < v)
    x[!delta] <- v[!delta]
    
    out <- data.frame(time = x, status = as.numeric(delta))
    dj <- n_each_interval(out$time[out$status == 1], p)
    if (all(dj > 2)) break
  }
  
  return (out)
}

# sum(logdata) in each I[j]
sum_logdata_pwexp <- function (x, p)
{
  k <- length(p)
  index <- index_each_interval(x, p)
  logdata <- vector(length = k)
  for (j in 1:k)
  {
    logdata[j] <- sum(x[index[[j]]] - p[j])
  }
  
  return (logdata)
}

# log_likelihood
loglik_pwexp <- function(theta, x, p)
{
  logC <- log(auxiliar_pwexp(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwexp(x, p)
  nj <- n_each_interval(x, p)
  l <- sum(nj * (log(theta) + logC)) - sum(theta * logdata)
  return (l)
}

# MLE
mle_pwexp <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x, p)
  nj <- n_each_interval(x, p)
  if (any(nj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(x[index[[j]]] - p[j])
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j]))
    alpha[j] <- nj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- nj[k] * sum((x[index[[k]]] - p[k])) ^ (-1)
  
  return(alpha)
}

# loglikelihood censored
loglik_cens_pwexp <- function(theta, x, p)
{
  logC <- log(auxiliar_pwexp(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwexp(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * log(theta) + nj * logC) - sum(theta * logdata)
  
  return (l)
}

# MLE to censored data
mle_cens_pwexp <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  if (any(dj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(x$time[index[[j]]] - p[j])
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j]))
    alpha[j] <- dj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- dj[k] * sum((x$time[index[[k]]] - p[k])) ^ (-1)
  
  return(alpha)
}

# profile_log_likelihood
profile_loglik_pwexp <- function(p, x)
{
  theta <- mle_pwexp(x, p)
  logC <- log(auxiliar_pwexp(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwexp(x, p)
  nj <- n_each_interval(x, p)
  l <- sum(nj * (log(theta) + logC)) - sum(theta * logdata)
  return (l)
}

# profile_log_cens_likelihood
profile_loglik_cens_pwexp <- function(p, x)
{
  theta <- mle_cens_pwexp(x, p)
  logC <- log(auxiliar_pwexp(p, theta))[-(length(p) + 1)]
  logdata <- sum_logdata_pwexp(x$time, p)
  nj <- n_each_interval(x$time, p)
  dj <- n_each_interval(x$time[x$status == 1], p)
  l <-
    sum(dj * log(theta) + nj * logC) - sum(theta * logdata)
  
  return (l)
}

diff_c <- function(x) c(x[1] - x[2], x[2] - x[3], x[3] - x[4])

# MLE BC
mle_pwexp_bc <- function (x, p)
{
  k <- length(p)
  alpha <- vector(length = k)
  index <- index_each_interval(x, p)
  nj <- n_each_interval(x, p)
  if (any(nj == 0)) {
    stop(paste("Error en la definicion de intervalos", "\n", ""))
  }
  
  for (j in 1:(k - 1))
  {
    a <- sum(x[index[[j]]] - p[j])
    b <- sum(nj[(j + 1):k] * (p[(j + 1)] - p[j]))
    alpha[j] <- nj[j] * (a + b) ^ (-1)
  }
  
  alpha[k] <- nj[k] * sum((x[index[[k]]] - p[k])) ^ (-1)
  
  #BC
  
  h_ij <- diag(-sum(nj)*diff_c(auxiliar_pwexp(p, alpha)) / alpha^2)
  h_ij1 <- diag(c(2*sum(nj) * diff_c(auxiliar_pwexp(p, alpha))[1] / alpha[1]^3, 0, 0))
  h_ij2 <- diag(c(0, 2*sum(nj)*diff_c(auxiliar_pwexp(p, alpha))[2] / alpha[2]^3, 0))
  h_ij3 <- diag(c(0, 0, 2*sum(nj)*diff_c(auxiliar_pwexp(p, alpha))[3] / alpha[3]^3))
  
  f_temp <- function(x) -sum(nj)*diff_c(auxiliar_pwexp(p, x)) / x^2
  temp <- jacobian(f_temp, alpha)
  h_ij_1 <- diag(temp[,1])
  h_ij_2 <- diag(temp[,2])
  h_ij_3 <- diag(temp[,3])
  
  A1 <- h_ij_1 - 0.5 * h_ij1
  A2 <- h_ij_2 - 0.5 * h_ij2
  A3 <- h_ij_3 - 0.5 * h_ij3
  A <- cbind(A1, A2, A3)
  
  Bias <- c(solve(h_ij) %*% A %*% fBasics::vec(solve(h_ij)))
  alpha <- alpha - Bias
  
  return(alpha)
}


