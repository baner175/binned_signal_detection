rm(list = ls())

library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)

real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)
eps <- 1e-3

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)

# SIGNAL DENSITY:
fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(dtrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd)*exp(x))
}
# SIGNAL CDF:
Fs <- function(x, mean = mean_sig, sd = sd_sig)
{
  return(ptrunc(exp(x), spec = 'norm', a = real_l, b = real_u,
                mean = mean, sd = sd))
}
# Finding the signal region around log(mean_sig) with mass 1-eps:
find_d <- function(d)
{
  pl <- Fs(log(mean_sig)-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(log(mean_sig)+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(log(mean_sig) - l,u - log(mean_sig)))
r <- sol$root
M_lower <- log(mean_sig) - r # lower bound of the signal region
M_upper <- log(mean_sig) + r # upper bound of the signal region

phys_data <- read.table('Fermi_LAT_physics.txt', header = TRUE)
x <- log(phys_data$x)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
ni <- hist(x, breaks = bin_ends, plot = FALSE)$count
N <- sum(ni)
############################ BASELINE MODEL ####################################
q_model <- function(alpha){
  q_mass <- (1/alpha)*((l+1)^(-alpha) - (u+1)^(-alpha))
  qi <- sapply(1:k, function(i){
    integrate(function(t){
      ((t+1)^(-alpha-1))/q_mass
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(ni*log(qi)))
}
alpha_hat <- nlminb(start = 0.01,
                    objective = q_model,
                    upper = 10, lower = 0)$par

q <- function(x){
  q_mass <- (1/alpha_hat)*((l+1)^(-alpha_hat) - (u+1)^(-alpha_hat))
  return(((x+1)^(-alpha_hat-1))/q_mass)
}
# constructing the proposal background g:
# means of the Gaussian components to mix with q_\alpha
mean1_in_g <- (M_lower + log(mean_sig))/2 
mean2_in_g <- (M_upper + log(mean_sig))/2
sig_fs <- sqrt(integrate(function(x) {(x^2)*fs(x)}, l, u)$value - integrate(function(x) {(x)*fs(x)}, l, u)$value^2)
sig_0 <- 3*sig_fs # SD for the gaussian components
# Defining proposal background g:
g <- function(x, lambda){
  phi_1 <- dtrunc(x, spec = 'norm', a = l, b = u,
                  mean = mean1_in_g, sd = sig_0)
  phi_2 <- dtrunc(x, spec = 'norm', a = l, b = u,
                  mean = mean2_in_g, sd = sig_0)
  return(
    lambda*(phi_1 + phi_2) + (1-2*lambda)*q(x)
  )
}

d_log_h <- function(t) -log(t+1)
E_q_d_log_h <- integrate(function(t) d_log_h(t)*q(t), l, u)$value
d_log_q <- function(t) d_log_h(t) - E_q_d_log_h
d2_log_q_int1 <- integrate(function(t) {
  (d_log_h(t)^2)*q(t)
}, l, u)$value
d2_log_q <- -d2_log_q_int1 + E_q_d_log_h^2
signal_search <- function(lambda){
  norm_S <- integrate(function(x){
    return(((fs(x)/g(x, lambda = lambda) - 1)^2)*g(x, lambda = lambda))
  }, l, u)$value |> sqrt()
  
  S0_xi <- sapply(xi, function(x){
    S_val <- (fs(x)/g(x, lambda = lambda) - 1)
    return(S_val/(norm_S^2))
  })
  theta0_hat <- sum(ni*S0_xi)/N
  d_log_q_xi <- sapply(xi, d_log_q)
  d_normS_sq <- -(1-2*lambda)*integrate(function(y){
    q <- q(y)
    fs <- fs(y)
    g <- g(y, lambda = lambda)
    d_log_q <- d_log_q(y)
    return(((fs/g)^2)*q*d_log_q)
  }, l, u)$value
  
  d_S0_xi <- sapply(xi, function(y){
    fs <- fs(y)
    q <- q(y)
    g <- g(y, lambda = lambda)
    d_log_q <- d_log_q(y)
    return(-((norm_S^2)*(fs/(g^2))*(1-2*lambda)*q*d_log_q + (fs/g-1)*d_normS_sq)/(norm_S^4))
  })
  V_hat <- sum((d_log_q_xi^2)*ni)/N
  J_hat <- -d2_log_q
  D_alh <- sum(d_S0_xi*ni)/N
  var_S0_F_hat <- sum((S0_xi^2)*ni)/N - theta0_hat^2
  sig_theta0_hat <- sqrt(
    var_S0_F_hat + (V_hat/(J_hat^2))*(D_alh^2) + 
      (2/J_hat)*D_alh*
      sum(S0_xi*d_log_q_xi*ni)/N
  )
  theta0_stat <- sqrt(N)*(theta0_hat-0)/sig_theta0_hat
  p_val <- pnorm(theta0_stat, lower.tail = FALSE)
  CI_95 <- theta0_hat + c(-1,1)*qnorm(0.975)*sig_theta0_hat/sqrt(N)
  return(c(theta0_hat,
           CI_95,
           p_val))
}

lambda_seq <- c(0.03, 0.05, 0.07)
res_sig_search <- sapply(lambda_seq, signal_search)
res_sig_search <- cbind(lambda_seq, t(res_sig_search))
res_sig_search <- data.frame(res_sig_search)
colnames(res_sig_search) <- c("lambda", "theta0beta_hat",
                              "CI_Lower", "CI_Upper", "p-value")
caption <- paste0("Signal Search Results")
kable(res_sig_search, format = "simple", digits = 10,
      booktabs = TRUE, escape = FALSE,
      caption = caption)
