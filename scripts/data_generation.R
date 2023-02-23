##############################

# Preparation-------------------------------------------------------------------
library(lme4)
library(lmerTest)
library(MASS)


# Design matrix ----------------------------------------------------------------
# n1: Size of the cluster.
# n2: Number of clusters.
# rho: Intraclass correlation.
# effect.size: Effect size
n1 <- 10
n2 <- 100
rho <- 0.11
effect_size <- 0.24
design_matrix <- expand.grid(n1, n2, rho, effect_size)

# Data generation---------------------------------------------------------------
# Set fixed and random effects
set.seed(26)
betas <- rnorm(3)
b_0 <- betas[1]
b_1 <- betas[2]
b_3 <- betas[3]
e_sd <- rnorm(1) # first level residual variance.
u0_sd <- rnorm(1) # random intercept.
u1_sd <- rnorm(1) # random slope.

cov <- rho * u0_sd * u1_sd
cox_matrix <- matrix(c(u0_sd^2, cov, cov, u1_sd^2), 2, 2)
u0 <- rep(rnorm(n2, sd = u0_sd))
e <- rnorm(n1 * n2, sd = e_sd)

# data frame
# cluster
cluster <- rep(1:n2, n1)
# treatment vs. control
sum_group <- 0
while (sum_group != n2 / 2) {
  group <- rbinom(n = n2, size = 1, prob = 0.5)
  sum_group <- sum(group)
}
