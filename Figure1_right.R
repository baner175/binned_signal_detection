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

dat <- read.table('Fermi_LAT_physics.txt', header = TRUE)$x
x <- log(dat)
N <- length(x)
k <- 1e2
bin_ends <- seq(l, u, length.out = k+1)
xi <- (bin_ends[1:k] + bin_ends[2:(k+1)])/2
ni <- ni <- sapply(1:k, function(i){
  sum((x > bin_ends[i])&(x <= bin_ends[i+1]))
})

qb_y_model <- function(beta){
  qb_mass <- (1/beta)*((l+1)^(-beta) - (u+1)^(-beta))
  qb_i <- sapply(x, function(t){
    ((t+1)^(-beta-1))/qb_mass
  })
  return(-sum(log(qb_i)))
}

beta_hat <- nlminb(start = 0.01,
                   objective = qb_y_model,
                   upper = Inf, lower = 0)$par

q <- function(x)
{
  qb_mass <- (1/beta_hat)*((l+1)^(-beta_hat) - (u+1)^(-beta_hat))
  ((x+1)^(-beta_hat-1))/qb_mass
}
# Figuring out (mu_s -d, mu_s + d) that covers an area of 1-eps:
find_d <- function(d)
{
  pl <- Fs(log(mean_sig)-d, mean = mean_sig, sd = sd_sig)
  pu <- Fs(log(mean_sig)+d, mean = mean_sig, sd = sd_sig)
  return(pu-pl-1+eps)
}

sol <- uniroot(find_d, lower = 0, upper = min(log(mean_sig) - l,u - log(mean_sig)))

r <- sol$root

M_lower <- log(mean_sig) - r
M_upper <- log(mean_sig) + r

# PROPOSAL BACKGROUND DENSITY:
mean1_in_gb <- (M_lower + log(mean_sig))/2
mean2_in_gb <- (M_upper + log(mean_sig))/2
# mean1_in_gb <- M_lower
# mean2_in_gb <- M_upper
sig_fs <- sqrt(integrate(function(x) {(x^2)*fs(x)}, l, u)$value - integrate(function(x) {(x)*fs(x)}, l, u)$value^2)
sd_in_gb <- 3*sig_fs

lambda_seq <- c(0, 0.01, 0.03, 0.05, 0.07)

######### Plotting Densities ###################################################

g <- function(y, lambda){
  q <- q(y)
  fs_val1 <- dtrunc(y, spec = 'norm', a = l, b = u,
                    mean = mean1_in_gb, sd = sd_in_gb)
  fs_val2 <- dtrunc(y, spec = 'norm', a = l, b = u,
                    mean = mean2_in_gb, sd = sd_in_gb)
  gb <- lambda*(fs_val1+fs_val2) + (1-2*lambda)*q
  return(gb)
}

op <- par(no.readonly = TRUE)
par(mar = c(5,5,3,2),
    mgp = c(2.5,1,0))
plot(x = xi, y = ni,
     main = '', 
     ylab = 'Counts',
     xlab = 'log(x)',
     col = 'grey30',
     pch = 16,
     cex.lab = 2,
     cex.axis = 2,
     ylim = c(0, max(ni)+10))

mycols <- c('cyan4', 'brown', 'darkgreen', 'orange', 'purple')
palette(mycols)
my_lty <- c(1,2,6,5,4)
for(j in 1:length(lambda_seq))
{
  curve(N*g(x, lambda = lambda_seq[j])*(u-l)/(k+1),
        l, u, add = TRUE, lwd = 4,
        col = mycols[j],
        lty = my_lty[j])
}
abline(v = c(M_lower, M_upper), col = 'grey', lty = 2, lwd = 4)
rect(
  xleft   = M_lower,
  xright  = M_upper,
  ybottom = -5,
  ytop    = max(ni)+15,
  col     = rgb(0, 1, 0, 0.2),
  border  = NA
)

legend(x = 0.95, y = 110,
       col = mycols,
       lty = my_lty, bty = 'n', lwd = 4,
       legend=c(TeX('$g_{\\hat{\\beta}}(\\lambda = 0) \\equiv q_{\\hat{\\alpha}}$'),
                TeX(sprintf('$g_{\\hat{\\beta}}(\\lambda = %f)$', lambda_seq[-1]))),
       cex = 2,
       x.intersp = 0.5,
       seg.len = 2.5,
       y.intersp = 1)
par(op)