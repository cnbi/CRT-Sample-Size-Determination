############# Sample Size Determination for Cluster Randomized Trials evaluating  #############
############################### informative hypotheses ######################################

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


SSD_crt_inform <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh, 
                         eta = 0.8, fixed = 'n2') {
    # Libraries
    library(lme4)
    library(bain)
    library(ggplot2)
    library(dplyr)
    
    #Functions
    source('data_generation.R')
    source('evaluation.R')
    
    # Starting values
    total.var <- 1
    var.u0 <- rho * total.var
    var.e <- total.var - var.u0
    iterations <- 1
    condition <- FALSE #condition fulfillment indicator 
    
    #Hypothesis
    hypothesis1 <- "Dintervention>Dcontrol"
    hypothesis2 <- "Dintervention<Dcontrol" 
    hypoth <- paste(hypothesis1, ";", hypothesis2)
    
    # Simulation
    while (condition == FALSE) {
        # If H1 is true
        data_crt <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, 
                                              eff.size, hypoth))
        print("Yay data ready!")
        
        # Evaluation of condition
        # Proportion
        prop.BF12 <- length(which(data_crt[, 'BF.12'] > BF.thresh)) / n.datasets 
        prop.BF21 <- length(which(data_crt[, 'BF.21'] < 1/BF.thresh)) / n.datasets # I am not sure of this, is it really necessary?
        # Evaluation
        ifelse(prop.BF12 > eta & prop.BF21 > eta, condition <- TRUE, condition <- FALSE)
        
        # If condition == FALSE then increase the sample.
        if (condition == FALSE) {
            print(c("Using cluster size: ", n1, " and number of clusters: ", n2, " eta:", prop.BF12))
            if (fixed == 'n1') {
                n2 = n2 + 2
            } else if (fixed == 'n2') {
                n1 = n1 + 1
            }
        }
        
        # Stop because is crazy the number of clusters, is not plausible
        iterations <- iterations + 1
        if (n2 == 1000) {
            break
        }
        
    }
    SSD_object <- list("n1" = n1,
                       "n2" = n2,
                       "Eta" = prop.BF12,
                       "data" = data_crt)
    # Output
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")
    cat("Hypotheses:", "\n")
    cat("H1:", hypothesis1, "\n")
    cat("H2:", hypothesis2, "\n")
    cat("Using cluster size = ", SSD_object$n1, " and number of clusters = ", SSD_object$n2, "\n")
    cat("P (BF.12 > BF.threshold | H.1) = ", SSD_object$Eta, "\n")
    return(SSD_object)
    
}


#TODO:
# - Check style with lintR
# - Think a better name for the condition object.
# - Add plots.This could be a function.
# - Add binary search algorithm.
# - Run a simulation to know see the performance under various conditions.
# - Error messages.(class of variables)


# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed are not TRUE both.

# Test -------------------------------------------------------------------------
start.time <- Sys.time()
a <- SSD_crt_inform(eff.size = 0.5, n.datasets = 100, rho = 0.1, BF.thresh = 3, fixed = 'n1')
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt_inform(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, fixed = "n1") #singularity
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


start.time <- Sys.time()
SSD_crt_inform(eff.size = 0.4, n.datasets = 15, rho = 0.01, BF.thresh = 6, fixed = "n1")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# This is weird, it seems like this needs less number of clusters.

start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.2, n.datasets = 20, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.2, n.datasets = 20, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.8, n.datasets = 20, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

