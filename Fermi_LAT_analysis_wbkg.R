rm(list = ls())

################################################################################
################################################################################
##                                                                            ##
##     Code for reproducing the analysis on the simulated data                ##
##     from Fermi-LAT in the presence of a background-only sample             ##
##                                                                            ##
################################################################################
################################################################################

####################### Loading the required packages ##########################
library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)
################################################################################

# Loading the data and transforming into log-scale
real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}

phys_data <- read.table('Fermi_LAT_physics.txt', header = TRUE)
bkg_data <- read.table('Fermi_LAT_bkg_only.txt', header = TRUE)
y <- log(bkg_data$x)
x <- log(phys_data$x)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
mi <- hist(y, breaks = bin_ends, plot = FALSE)$count
ni <- hist(x, breaks = bin_ends, plot = FALSE)$count
N <- sum(ni)
M <- sum(mi)

################################################################################
######## PARAMETRIC MODEL: TRUNCATED EXPONENTIAL ###############################
################################################################################

g_model <- function(beta){ # defining log-likelihood based g_beta
  gi <- sapply(1:k, function(i){
    integrate(function(x){
      dtrunc(x, spec = 'exp', a = l, b = u,
             rate = beta)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(mi*log(gi)))
}
beta_hat <- nlminb(start = 0.01,
                   objective = g_model,
                   upper = 10, lower = 0)$par
# Defining proposal background g at MLE beta_hat:
g <- function(x) dtrunc(x, spec = 'exp', rate = beta_hat, a = l, b = u)

norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g) # ||S||_{G_\beta} at \beta = \hat{\beta}
}, l, u)$value |> sqrt()
S0_vec <- sapply(xi, function(t){ # Evaluating S_0 on the bin-mids
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- sum(S0_vec*ni)/N; delta0_hat <- sum(S0_vec*mi)/M
eta_hat_exp <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta
test_num <- sqrt(M*N)*eta_hat_exp # numerator of the test statistic

###### Now we shall compute the components involved ###############
###### in the denominator of the test statistic ###################

# derivative of (||S||_G)^2 w.r.t. \beta
d_normS_sq <- -integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- (1/beta_hat) - t - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  return((fs^2)*d_log_g/g)
}, l, u)$value
d_S0 <- function(t){ # derivative of S_0 w.r.t. \beta
  fs <- fs(t)
  g <- g(t)
  d_log_g <- (1/beta_hat) - t - 
    (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u))
  
  return(-((norm_S^2)*(fs/g)*d_log_g + (fs/g-1)*d_normS_sq)/(norm_S^4))
}
d_log_g_xi <- sapply(xi, function(t){ # evaluating d_\beta log(g) on the bin-mids
  return(
    (1/beta_hat) - t - 
      (u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)) 
  )
})
# double derivative of log g_\beta w.r.t. \beta:
d2_log_g<- -(1/(beta_hat^2)) - (((l^2)*exp(-beta_hat*l) - (u^2)*exp(-beta_hat*u))/(exp(-beta_hat*l) - exp(-beta_hat*u)) - 
                                  ((u*exp(-beta_hat*u) - l*exp(-beta_hat*l))/(exp(-beta_hat*l) - exp(-beta_hat*u)))^2)

d_theta0 <- sum(sapply(xi, d_S0)*ni)/N # derivative of \hat\theta_0 w.r.t. \beta
d_delta0 <- sum(sapply(xi, d_S0)*mi)/M # derivative of \hat\delta_0 w.r.t. \beta
d_theta_T <- 1/(1-delta0_hat) # derivative of T w.r.t. first component
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2) # derivative of T w.r.t. second component
# components for the denominator
cov_term <- sum(d_log_g_xi*S0_vec*mi)/M
V_hat <- sum((d_log_g_xi^2)*mi)/M
J_hat <- -sum(mi*d2_log_g)/M
var_S0_F_hat <- sum((S0_vec^2)*ni)/N - (theta0_hat^2)
var_S0_Fb_hat <- sum((S0_vec^2)*mi)/M - (delta0_hat^2)
denom1 <- M*var_S0_F_hat*(d_theta_T^2)
denom2 <- N*var_S0_Fb_hat*(d_delta_T^2)
denom3 <- N*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*N*(1/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)
test_denom <- sqrt(denom1 + denom2 + denom3 + denom4) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_exp <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_exp <- test_denom/(sqrt(M+N))
std_err_exp <- sig_hat_exp*sqrt((M+N)/(M*N)) # standard error
ci_95_exp <- eta_hat_exp + c(-1,1)*qnorm(0.975)*std_err_exp # 95% CI

################################################################################
######## PARAMETRIC MODEL: TAIL OF TRUNCATED GAUSSIAN  #########################
################################################################################

mu_in_g <- -1; sigma_factor_in_g <- 2# so that g_\beta(x) \propto exp(-((x+1)^2)/(4*\beta))
# below we replicate the same code as above with necessary adjustments for the new background model
g_model <- function(beta){
  gi <- sapply(1:k, function(i){
    integrate(function(x){
      dtrunc(x, spec = 'norm',
             mean = mu_in_g,
             sd = sqrt(sigma_factor_in_g*beta),
             a = l, b = u)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(mi*log(gi)))
}
beta_hat <- nlminb(start = 0.01,
                   objective = g_model,
                   upper = 10, lower = 0)$par
g <- function(t) dtrunc(t, spec = 'norm',
                        mean = mu_in_g,
                        sd = sqrt(sigma_factor_in_g*beta_hat),
                        a = l, b = u)
d_log_h <- function(t) ((t-mu_in_g)^2)/(2*sigma_factor_in_g*(beta_hat^2))
d2_log_h <- function(t) (-(t-mu_in_g)^2)/(sigma_factor_in_g*(beta_hat^3))
E_g_d_log_h <- integrate(function(t) d_log_h(t)*g(t), l, u)$value
d2_log_g_int_1 <- integrate(function(y){
  d2_log_h <- d2_log_h(y)
  d_log_h <- d_log_h(y)
  d2_h_by_h <- d2_log_h + d_log_h^2
  g <- g(y)
  return(d2_h_by_h*g)
}, l, u)$value
d2_log_g_int_2 <- integrate(function(y){
  d_log_h <- d_log_h(y)
  g <- g(y)
  return(d_log_h*g)
}, l, u)$value
d_log_g <- function(t) d_log_h(t) - E_g_d_log_h
d2_log_g <- function(t){
  d2_log_h <- d2_log_h(t)
  int_val <- d2_log_g_int_1 - (d2_log_g_int_2^2)
  return(d2_log_h - int_val)
}
norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt()
S0_vec <- sapply(xi, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- sum(S0_vec*ni)/N; delta0_hat <- sum(S0_vec*mi)/M
eta_hat_GT <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta 
test_num <- sqrt(M*N)*eta_hat_GT # numerator of the test statistic
d_normS_sq <- -integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- d_log_g(t)
  return((fs^2)*d_log_g/g)
}, l, u)$value

d_S0 <- function(t){
  fs <- fs(t)
  g <- g(t)
  d_log_g <- d_log_g(t)
  return(-((norm_S^2)*(fs/g)*d_log_g + (fs/g-1)*d_normS_sq)/(norm_S^4))
}
d_log_g_xi <- sapply(xi, d_log_g)
d2_log_g_xi <- sapply(xi, d2_log_g)
d_theta0 <- sum(ni* sapply(xi, d_S0))/N
d_delta0 <-  sum(mi* sapply(xi, d_S0))/M
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
# calculating the denominator of the test statistic
cov_term <- sum(mi*d_log_g_xi*S0_vec)/M
V_hat <- sum((d_log_g_xi^2)*mi)/M
J_hat <- -sum(d2_log_g_xi*mi)/M
var_S0_F_hat <- sum((S0_vec^2)*ni)/N - (theta0_hat^2)
var_S0_Fb_hat <- sum((S0_vec^2)*mi)/M - (delta0_hat^2)
denom1 <- M*var_S0_F_hat*(d_theta_T^2)
denom2 <- N*var_S0_Fb_hat*(d_delta_T^2)
denom3 <- N*(V_hat/(J_hat^2))*((d_theta_T*d_theta0 + d_delta_T*d_delta0)^2)
denom4 <- 2*N*(1/J_hat)*d_delta_T*cov_term*(d_theta_T*d_theta0 + d_delta_T*d_delta0)
test_denom <- sqrt(denom1 + denom2 + denom3 + denom4) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_GT <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_GT <- test_denom/(sqrt(M+N))
std_err_GT <- sig_hat_GT*sqrt((M+N)/(M*N)) # std error
ci_95_GT <- eta_hat_GT + c(-1,1)*qnorm(0.975)*std_err_GT # p-value

################################################################################
##############  UNIFORM BACKGROUND (NO PARAMETERS INVOLVED)   ##################
################################################################################

g <- function(t) dunif(t, min = l, max = u) # uniform proposal background
norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt()
# similar code as above for the uniform background
S0_vec <- sapply(xi, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- sum(ni*S0_vec)/N; delta0_hat <- sum(mi*S0_vec)/M
eta_hat_unif <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta
test_num <- sqrt(M*N)*eta_hat_unif # numerator of the test statistic
# calculating the denominator of the test statistic
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
var_S0_F_hat <- sum((S0_vec^2)*ni)/N - (theta0_hat^2)
var_S0_Fb_hat <- sum((S0_vec^2)*mi)/M - (delta0_hat^2)
denom1 <- M*var_S0_F_hat*(d_theta_T^2)
denom2 <- N*var_S0_Fb_hat*(d_delta_T^2)
test_denom <- sqrt(denom1 + denom2) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_unif <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_unif <- test_denom/(sqrt(M+N))
std_err_unif <- sig_hat_unif*sqrt((M+N)/(M*N)) # std error
ci_95_unif <- eta_hat_unif + c(-1,1)*qnorm(0.975)*std_err_unif # 95% CI

################################################################################
#####################  EXP BACKGROUND FIXED RATE  ##############################
################################################################################

g <- function(t) dtrunc(t, spec = 'exp', a = l, b = u,
                        rate = 0.5) # proposal background
norm_S <- integrate(function(t){
  fs <- fs(t)
  g <- g(t)
  return(((fs/g - 1)^2)*g)
}, l, u)$value |> sqrt()
# similar code as above for the relevant exponential background
S0_vec <- sapply(xi, function(t){
  fs <- fs(t)
  g <- g(t)
  S_val <- (fs/g - 1)
  return(S_val/(norm_S^2))
})
theta0_hat <- sum(ni*S0_vec)/N; delta0_hat <- sum(mi*S0_vec)/M
eta_hat_expf <- (theta0_hat-delta0_hat)/(1-delta0_hat) # estimate of eta
test_num <- sqrt(M*N)*eta_hat_expf # numerator of the test statistic
# calculating the denominator of the test statistic
d_theta_T <- 1/(1-delta0_hat)
d_delta_T <- (theta0_hat - 1)/((1-delta0_hat)^2)
var_S0_F_hat <- sum((S0_vec^2)*ni)/N - (theta0_hat^2)
var_S0_Fb_hat <- sum((S0_vec^2)*mi)/M - (delta0_hat^2)
denom1 <- M*var_S0_F_hat*(d_theta_T^2)
denom2 <- N*var_S0_Fb_hat*(d_delta_T^2)
test_denom <- sqrt(denom1 + denom2) # denominator of the test statistic
test_stat <- test_num/test_denom # test statistic
p_val_expf <- pnorm(test_stat, lower.tail = FALSE) # p-value
sig_hat_expf <- test_denom/(sqrt(M+N))
std_err_expf <- sig_hat_expf*sqrt((M+N)/(M*N)) # std error
ci_95_expf <- eta_hat_expf + c(-1,1)*qnorm(0.975)*std_err_expf # 95% CI

########### Table of results ###################################################
df_res <- data.frame('Proposal_bkg' = c('Exp(β)', 'Gaussian-tail(β)', 'Exp(0.5)', 'U[l,u]'),
                     'eta_hat' = c(eta_hat_exp, eta_hat_GT, eta_hat_expf, eta_hat_unif),
                     'CI_Lower' = c(ci_95_exp[1], ci_95_GT[1], ci_95_expf[1], ci_95_unif[1]),
                     'CI_Upper' = c(ci_95_exp[2], ci_95_GT[2], ci_95_expf[2], ci_95_unif[2]),
                     'p_val' = c(p_val_exp, p_val_GT, p_val_expf, p_val_unif))
caption <- paste0(sprintf("Signal Search Results (# of bins = %d)",k))
kable(df_res, format = "simple", digits = 10,
      booktabs = TRUE, escape = FALSE,
      caption = caption)

