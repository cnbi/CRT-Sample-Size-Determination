############# Sample Size Determination for Cluster Randomized Trials #############

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description 
#'@arguments
## eff.size: Effect size
## n.datasets: Number of datasets that the user wants to generate to determine the sample size.
## var.u0: Between cluster variance.
## BFthresh: Threshold to choose.
## eta: Posterior probability that is expected. Related to the Bayesian error probability.
## plots: Logical argument. If TRUE the function print the plost.
## hypothesis: The hypothesis that is going to be tested.


SSD_crt <- function(eff.size = , n.datasets = 100, var.u0 = , BFthresh = , eta = , plots = TRUE, hypothesis = ) {
    # Libraries
    library(lme4)
    library(bain)
    library(ggplot2)
    library(dplyr)
    library(scales)
    
    #Functions
    source('data_generation.R')
    
    # Starting values
    n1 <- 15
    n2 <- 30
    total.var <- 1
    mean.ctrl <- 0
    rho <- var.u0/total.var #
    var.e <- total.var - var.u0
    
    # while starts here
    # If H1 is true
    results.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho, eff.size))
    
    # If H0 is true
    eff.size <- 0
    results.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho, eff.size))
    
    
    
    
    
    #Evaluation of condition
    
    # output
    
    
}