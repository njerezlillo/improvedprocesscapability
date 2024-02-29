library(ggQC)
library(dplyr)
library(maxLik)
library(xtable)
library(survival)
library(numDeriv)
library(ggfortify)
source("capability.R")
source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwexp.R")

df_2 <- read.csv("Application/dataset2.csv")[, 2]
length(df_2)
fit <- survfit(Surv(df_2) ~ 1)
tmin <- min(df_2)
n <- length(df_2)
shapiro.test(df_2)

# Change points estimation ------------------------------------------------

target <- function(y) profile_loglik_pwexp_bc(p = c(tmin, y), df_2)

# k = 1

A_1 <- matrix(1, ncol = 1, byrow = T)
d_1 <- c(-294)
fit_p_1 <- maxSANN(target, start = 300,
                   constraints = list(ineqA = A_1, ineqB = d_1))
fit_p_1 <- c(tmin, fit_p_1$estimate)
fit_alpha_1 <- mle_pwexp_bc(df_2, fit_p_1)

# k = 2

A_2 <- matrix(c(1, 0, -1, 1, 0, -1), ncol = 2, byrow = T)
d_2 <- c(-294, -1, 320)
fit_p_2 <- maxSANN(target, start = c(300, 305),
                   constraints = list(ineqA = A_2, ineqB = d_2))
fit_p_2 <- c(tmin, fit_p_2$estimate)
fit_alpha_2 <- mle_pwexp_bc(df_2, fit_p_2)

# k = 3

A_3 <- matrix(c(1, 0, 0, -1, 1, 0, 0, -1, 1, 0, 0, -1), ncol = 3, byrow = T)
d_3 <- c(-294, -1, -1, 320)
fit_p_3 <- maxSANN(target, start = c(300, 303, 306),
                   constraints = list(ineqA = A_3, ineqB = d_3))
fit_p_3 <- c(tmin, fit_p_3$estimate)
fit_alpha_3 <- mle_pwexp_bc(df_2, fit_p_3)

# AIC ---------------------------------------------------------------------

AIC_1 <- 2 * 3 - 2 * loglik_pwexp(fit_alpha_1, df_2, fit_p_1)
AIC_2 <- 2 * 5 - 2 * loglik_pwexp(fit_alpha_2, df_2, fit_p_2)
AIC_3 <- 2 * 7 - 2 * loglik_pwexp(fit_alpha_3, df_2, fit_p_3)

BIC_1 <- log(n) * 3 - 2 * loglik_pwexp(fit_alpha_1, df_2, fit_p_1)
BIC_2 <- log(n) * 5 - 2 * loglik_pwexp(fit_alpha_2, df_2, fit_p_2)
BIC_3 <- log(n) * 7 - 2 * loglik_pwexp(fit_alpha_3, df_2, fit_p_3)

gof <- data.frame(
  "k" = 1:3,
  AIC = round(c(AIC_1, AIC_2, AIC_3), 2),
  BIC = round(c(BIC_1, BIC_2, BIC_3), 2)
)

gof

# print(xtable(gof), include.rownames = F)

# Kolmogrov-Smirnov Test --------------------------------------------------

F_exp <- function(x) {
  1 - spwexp(x, fit_p_1, fit_alpha_1)
}

F_exp <- Vectorize(F_exp)

goftest::cvm.test(df_2, "F_exp")

F_exp <- function(x) {
  1 - spwexp(x, fit_p_2, fit_alpha_2)
}

F_exp <- Vectorize(F_exp)

goftest::cvm.test(df_2, "F_exp")

F_exp <- function(x) {
  1 - spwexp(x, fit_p_3, fit_alpha_3)
}

F_exp <- Vectorize(F_exp)

goftest::cvm.test(df_2, "F_exp")

# Estimates ---------------------------------------------------------------

se <- fit_alpha_1 / sqrt(n * diff_c(auxiliar_pwexp(fit_p_1, fit_alpha_1)))
ic_l <- round(fit_alpha_1 - 1.96 * se, 2)
ic_u <- round(fit_alpha_1 + 1.96 * se, 2)

est <- data.frame(j = 1:2,
                  "estimation" = round(fit_alpha_1, 2),
                  "CI" = paste0("(", ic_l, ", ", ic_u, ")"))
est

# print(xtable(est), include.rownames = F)

# Cpk ---------------------------------------------------------------------

M <- qpwexp(0.5, fit_p_1, fit_alpha_1)
Lp <- qpwexp(0.00135, fit_p_1, fit_alpha_1)
Up <- qpwexp(0.99865, fit_p_1, fit_alpha_1)
LSL <- 285
USL <- 315
H <- diag(se^2)
Cpk_real <- min((USL - M)/(Up - M), (M - LSL)/(M - Lp))
EE_Cpk <- md_Cpk(fit_alpha_1, H, fit_p_1, LSL, USL)

Cpk_real
Cpk_real - 1.96 * sqrt(EE_Cpk)
Cpk_real + 1.96 * sqrt(EE_Cpk)

# Cpm ---------------------------------------------------------------------

Tc <- 300
Cpm_real <- (USL - LSL)/(3 * sqrt(((Up - Lp)/3)^2 + (M - Tc)^2)) 
EE_Cpm <- md_Cpm(fit_alpha_1, H, fit_p_1, LSL, USL, Tc)

Cpm_real
Cpm_real - 1.96 * sqrt(EE_Cpm)
Cpm_real + 1.96 * sqrt(EE_Cpm)

# Cpm* --------------------------------------------------------------------

temp <- min(USL - Tc, Tc - LSL)
CpmA_real <- temp/(3 * sqrt(((Up - Lp)/3)^2 + (M - Tc)^2)) 
EE_CpmA <- md_CpmA(fit_alpha_1, H, fit_p_1, LSL, USL, Tc)

CpmA_real
CpmA_real - 1.96 * sqrt(EE_CpmA)
CpmA_real + 1.96 * sqrt(EE_CpmA)

# Cpmk --------------------------------------------------------------------

temp1 <- (USL - M)/(3 * sqrt(((Up - M)/3)^2 + (M - Tc)^2)) 
temp2 <- (M - LSL)/(3 * sqrt(((M - Lp)/3)^2 + (M - Tc)^2))
Cpmk_real <- min(temp1, temp2)
EE_Cpmk <- md_Cpmk(fit_alpha_1, H, fit_p_1, LSL, USL, Tc) 

Cpmk_real
Cpmk_real - 1.96 * sqrt(EE_Cpmk)
Cpmk_real + 1.96 * sqrt(EE_Cpmk)

# Figures -----------------------------------------------------------------

### Histogram

histogram <- ggplot(data.frame(df_2), aes(x = df_2)) +
  geom_histogram(
    binwidth = 3.5,
    color = "black",
    fill = "lightblue",
    boundary = 60
  ) +
  geom_vline(
    xintercept = c(285 - 0.5, 300 - 0.5, 315 - 0.5),
    color = c("red", "darkgreen", "red"),
    linetype = "dashed",
    size = 0.8
  ) +
  scale_x_continuous(
    breaks = c(285 - 0.5, 300 - 0.5, 315 - 0.5),
    labels = c("285", "300", "315")
  ) +
  scale_y_continuous(breaks = seq(0, 50, by = 5)) +
  labs(x = "Capacitance (μF)", y = "Frequency") +
  theme_bw() +
  ggtitle("a)") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15),
        plot.title = element_text (size = 15))

### KM

S_exp <- function(x) {
  spwexp(x, fit_p_1, fit_alpha_1)
}

S_exp <- Vectorize(S_exp)

p_df <- data.frame(p = fit_p_1,
                   e = S_exp(fit_p_1))

reliability <- 
  autoplot(fit, conf.int.fill = 'dodgerblue3', conf.int.alpha = 0.1) + theme_bw() +
  stat_function(fun=S_exp, col = "dodgerblue3", size = 0.7) + 
  geom_point(aes(x = p, y = e), col = "dodgerblue3", data = p_df, size = 1.2) +
  labs(x = "time", y = "Reliability") + xlim(c(290, 320)) +
  ggtitle("b)") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15),
        plot.title = element_text (size = 15))

### Cox-Snell

ei <- -log(S_exp(df_2)) 
km_ei <- survfit(Surv(ei) ~ 1)

residuals <-
  ggplot() + aes(x = km_ei$surv, y = exp(-km_ei$time)) + 
  geom_point() + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, color = "red", size = 1.2) +
  labs(x = "S(ei): Kaplan-Meier", y = "S(ei):standard exponential") +
  ggtitle("c)") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15),
        plot.title = element_text (size = 15))

### Sensibility

GCD <- vector(length = n)
Angle <- vector(length = n)

for (i in 1:n) {
  alpha_aux <- mle_pwexp(df_2[-i], fit_p_1)
  theta_i <- (alpha_aux - fit_alpha_1)
  GCD[i] <- t(theta_i) %*% solve(H) %*% theta_i
  Angle[i] <- (atan(alpha_aux[2] / alpha_aux[1]) -
               atan(fit_alpha_1[2] / fit_alpha_1[1])) * (180 / pi)
}

gcd_df <- data.frame(x = 1:n, y = GCD) %>% arrange(desc(y))

sensibility <- ggplot(gcd_df) + 
  aes(x, y) + geom_point(pch = 1) + 
  geom_text(data = gcd_df[1,], 
            aes(x = x, y = y + 0.02, label = x[1]),
            size = 5) +
  labs(x = "Index", y = "Generalized Cook Distance") + 
  theme_bw() +
  ggtitle("d)") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text (size = 15),
        plot.title = element_text (size = 15))

Angle[100]

g_km <- gridExtra::grid.arrange(histogram, reliability, residuals, sensibility, ncol = 2)
ggsave(plot = g_km, "./Application/Fig2.pdf", height = 12, width = 15)
