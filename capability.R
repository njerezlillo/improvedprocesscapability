
#rpwexp <- function(n, p, alpha)
#{
#  q <- Vectorize(function(x) qpwexp(x, p, alpha), "x")
#  while(TRUE)
#  {
#    x <- q(runif(n))
#    nj <- n_each_interval(x, p)
#    if (any(nj < 1)) {
#      x[1] <- p[2] - 0.1
#      break
#    }
#  }
#  
#  return (x)
#}

md_Cpk <- function (mod, cov, p, LSL, USL) {
  myfun <- function (z) {
    Up_aux <- qpwexp(0.99865, p, z)
    Lp_aux <- qpwexp(0.00135, p, z)
    M_aux <- qpwexp(0.5, p, z)
    min((USL - M_aux)/(Up_aux - M_aux),
        (M_aux - LSL)/(M_aux - Lp_aux))
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}

md_Cpmk <- function (mod, cov, p, LSL, USL, Tc) {
  myfun <- function (z) {
    Up_aux <- qpwexp(0.99865, p, z)
    Lp_aux <- qpwexp(0.00135, p, z)
    M_aux <- qpwexp(0.5, p, z)
    aux1 <- (USL - M_aux)/(3 * sqrt(((Up_aux - M_aux)/3)^2 + (M_aux - Tc)^2))
    aux2 <- (M_aux - LSL)/(3 * sqrt(((M_aux - Lp_aux)/3)^2 + (M_aux - Tc)^2))
    min(aux1, aux2)
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}

md_Spmk <- function (mod, cov, p, LSL, USL, Tc) {
  myfun <- function (z) {
    mu_aux <- 
      integrate(Vectorize(function(x) x * dpwexp(x, p, z)), 0, Inf)$value
    var_aux <- 
      integrate(Vectorize(function(x) (x - mu_aux)^2 * dpwexp(x, p, z)), 0, Inf)$value
    
    qnorm((1 + ppwexp(LSL, p, z) - ppwexp(USL, p, z))/2) /
      (3 * sqrt(1 + (mu_aux - Tc)^2/var_aux))
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}

md_Cpm <- function (mod, cov, p, LSL, USL, Tc) {
  myfun <- function (z) {
    Up_aux <- qpwexp(0.99865, p, z)
    Lp_aux <- qpwexp(0.00135, p, z)
    M_aux <- qpwexp(0.5, p, z)
    (USL - LSL)/(3 * sqrt(((Up_aux - Lp_aux)/3)^2 + (M_aux - Tc)^2))
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}

md_CpmA <- function (mod, cov, p, LSL, USL, Tc) {
  myfun <- function (z) {
    Up_aux <- qpwexp(0.99865, p, z)
    Lp_aux <- qpwexp(0.00135, p, z)
    temp <- min(USL - Tc, Tc - LSL)
    M_aux <- qpwexp(0.5, p, z)
    temp/(3 * sqrt(((Up_aux - Lp_aux)/3)^2 + (M_aux - Tc)^2))
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}
