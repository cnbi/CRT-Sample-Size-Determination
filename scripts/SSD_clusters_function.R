############# Sample Size Determination for Cluster Randomized Trials #############

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description 
#'@arguments
## eff.size: Effect size
## n.datasets: Number of data sets that the user wants to generate to determine the sample size.
## var.u0: Between cluster variance.
## BF.thresh: Threshold to choose.
## eta: Posterior probability that is expected. 1-eta = Bayesian error probability.
## plots: Logical argument. If TRUE the function print the plots.
## hypothesis: The hypothesis that is going to be tested.


SSD_crt <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 100, rho, BF.thresh, 
                    eta = 0.8, plots = TRUE, hypothesis, n1.fixed = TRUE,
                    n2.fixed = TRUE) {
    # Libraries
    library(lme4)
    library(bain)
    library(ggplot2)
    library(dplyr)
    library(scales)
    
    #Functions
    source('data_generation.R')
    source('extract_results.R')
    source('evaluation.R')
    
    # Starting values
    total.var <- 1
    mean.ctrl <- 0
    var.u0 <- rho * total.var
    var.e <- total.var - var.u0
    
    # Names of table with results
    names.table <- c('Dcontrol', 'Dtreatment', 'BF.12', 'BF.21', 'BFc.1', 'BFc.2', 'PMP.1',
                     'PMP.2', 'PMPc.1', 'PMPc.2', 'PMPc.c')
    condition <- FALSE #condition fulfillment indicator
    
    while (condition == FALSE) {
        # If H1 is true
        data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho, eff.size))
        results.H1 <- do.call(extract_results, list(n.datasets, data.H1, names.table = names.table))
        
        # If H0 is true
        eff.size <- 0
        data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho, eff.size))
        results.H0 <-  do.call(extract_results, list(n.datasets, data.H0, names.table = names.table))
        
        #Evaluation of condition
        condition <- eval_thresh(results.H0 = results.H0, results.H1 = results.H1, 
                                 BF.thresh = BF.thresh, n.datasets = n.datasets, condition = condition, eta = eta)
        
        # If condition =FALSE then increase the sample. Which n increase?
        if (condition == FALSE){
            if (n1.fixed == TRUE) {
                n2 = n2 + 2
            } else (n2.fixed == TRUE){
                n1 = n1 + 1
            }
        }
        
    }
    
    # output
    #medians.results <- apply(results, 2, median) #returns a vector for plots
    
    
}


#TODO:
    # - Add ICC calculation
    # - Think: If n1.fixed and n2.fixed are FALSE, then what?
    # - Think: How is going to be the output?
    # - Think: Should I change BF12 to BF01 and BF21 to BF10?
    # - Check style with lintR
    # - Talk w/ Uli about the name of the arguments → They should be the same.
    # - Test evaluation function.
    # - Think a better name for the condition.
    # - Hypothesis should be entered by researcher. Change this in data_generation function.
    # - Add plots.This could be a function.
    # - Add binary search algorithm.

# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed aro not TRUE both.