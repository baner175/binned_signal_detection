rm(list = ls())

library(VGAM)
library(truncdist)
library(latex2exp)
library(knitr)
library(kableExtra)
library(ggplot2)

real_l <- 1; real_u <- 35
l <- log(real_l); u <- log(real_u)

mean_sig <- 3.5; sd_sig <- sqrt(0.01*3.5^2)
eps <- 1e-3
mu_in_g <- -1; sigma_factor_in_g <- 2

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

# True background:
fb <- function(x){
  return(dtrunc(x, spec = 'exp', a = l, b = u,
                rate = 1.4))
}

phys_data <- read.table('Fermi_LAT_physics.txt', header = TRUE)
bkg_data <- read.table('Fermi_LAT_bkg_only.txt', header = TRUE)
y <- log(bkg_data$x)
x <- log(phys_data$x)

N <- length(x); M <- length(y)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2

ni <- sapply(1:k, function(i){
  sum((x > bin_ends[i])&(x <= bin_ends[i+1]))
})
mi <- sapply(1:k, function(i){
  sum((y > bin_ends[i])&(y <= bin_ends[i+1]))
})

g_GT_model <- function(beta){
  g_i <- sapply(1:k, function(i){
    integrate(function(t){
      dtrunc(t, spec = 'norm', mean = mu_in_g,
             sd = sqrt(sigma_factor_in_g*beta),
             a = l, b = u)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(mi*log(g_i)))
}

beta_GT_hat <- nlminb(start = 0.01,
                   objective = g_GT_model,
                   upper = 10, lower = 0)$par

g_GT <- function(t) dtrunc(t, spec = 'norm',
                         mean = mu_in_g,
                         sd = sqrt(sigma_factor_in_g*beta_GT_hat),
                         a = l, b = u)

g_Exp_model <- function(beta){
  g_i <- sapply(1:k, function(i){
    integrate(function(t){
      dtrunc(t, spec = 'exp', rate = beta,
             a = l, b = u)
    }, bin_ends[i], bin_ends[i+1])$value
  })
  return(-sum(mi*log(g_i)))
}

beta_Exp_hat <- nlminb(start = 0.01,
                      objective = g_Exp_model,
                      upper = Inf, lower = 0)$par

g_Exp <- function(t) dtrunc(t, spec = 'exp', rate = beta_Exp_hat,
                             a = l, b = u)

g_Expf <- function(t) dtrunc(t, spec = 'exp', rate = 0.5,
                            a = l, b = u)



op <- par(no.readonly = TRUE)
par(mar = c(5,5,3,2),
    mgp = c(2.5,1,0))
plot(y = mi, x = xi,
     pch = 16, col = 'grey30',
     ylim = c(0, max(mi)+10),
     xlab = 'log(x)',
     ylab = 'Counts',
     cex.lab = 2,
     cex.axis = 2,)
curve(M*fb(x)*(u-l)/(k+1), col = 'black',
      add = TRUE,
      lty = 1, lwd = 4)
curve(M*g_Exp(x)*(u-l)/(k+1), col = 'blue', 
      add = TRUE,
      lty = 5, lwd = 4)
curve(M*g_GT(x)*(u-l)/(k+1), col = 'brown',
      add = TRUE,
      lty = 2, lwd = 4)
curve(M*g_Expf(x)*(u-l)/(k+1), col = 'purple', 
      add = TRUE,
      lty = 4, lwd = 4)
curve(M*dunif(x, l, u)*(u-l)/(k+1), col = 'red', 
      add = TRUE,
      lty = 6, lwd = 4)

legend(x = 0.85, y = 225,
       legend = c(
         TeX('$f_b(x)$'),
         TeX('$g_{\\hat{\\beta}}(x) \\propto \\exp(-\\hat{\\beta}x)$'),
         TeX('$g_{\\hat{\\beta}}(x) \\propto \\exp\\left{-\\frac{(x+1)^2}{2\\hat{\\beta}}\\right}$'),
         TeX('$g(x) \\propto \\exp(-0.5x)$'),
         TeX('$g(x) \\propto 1$')
       ),
       bty = 'n',
       cex = 1.7,
       lty = c(1,5,2,4,6),
       lwd = 4,
       col = c('black', 'blue', 'brown', 'purple', 'red'),
       y.intersp = 0.5,
       seg.len = 2.5
)
par(op)


