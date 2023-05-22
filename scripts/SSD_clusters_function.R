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


SSD_crt_null <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh, 
                    eta = 0.8, fixed = 'n2') {
    # Libraries
    library(lme4)
    library(bain)
    library(ggplot2)
    library(dplyr)
    
    #Functions
    source('data_generation.R')
    # source('print.null')
    
    # Starting values
    total.var <- 1
    var.u0 <- rho * total.var
    var.e <- total.var - var.u0
    iterations <- 1
    eff.size0 <- 0
    condition <- FALSE #condition fulfillment indicator
    
    #Hypotheses
    hypothesis1 <- "Dintervention>Dcontrol"
    null <- "Dintervention=Dcontrol"
    hypoth <- paste(null, ";", hypothesis1)
    
    # Simulation and evaluation of condition 
    while (condition == FALSE) {
        # If H1 is true
        data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, 
                                              eff.size))
        colnames(data.H1) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
        
        print("Yay data ready!")
        # If H0 is true
        data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e,
                                              eff.size = eff.size0))
        colnames(data.H0) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
        
        #Evaluation of condition
        # Proportion
        prop.BF01 <- length(which(data.H0[, 'BF.01'] > BF.thresh)) / n.datasets 
        prop.BF10 <- length(which(data.H1[, 'BF.10'] > BF.thresh)) / n.datasets 
        # Evaluation
        ifelse(prop.BF01 > eta & prop.BF10 > eta, condition <- TRUE, condition <- FALSE)
        #browser()
        #Output
        print(c())
        print("Yay evaluation works!")
        # If condition =FALSE then increase the sample. Which n increase?
        if (condition == FALSE) {
            print("Using cluster size:", n1, "and number of clusters:", n2, 
                  "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10, sep = " ")
            if (fixed == 'n1') {
                n2 = n2 + 2
            } else if (fixed == 'n2') {
                n1 = n1 + 1
            }
        }
        iterations <- iterations + 1
        print(c("n1: ", n1, " n2: ", n2))
        if (n2 == 1000) {
            break
        }
        
    }
    
}

# output
print(c("Cluster size: ", n1, "Number of clusters: ", n2))




#TODO:
    # - Think: If n1.fixed and n2.fixed are FALSE, then what?
    # - Think: How is going to be the output?
    # - Think: Should I change BF12 to BF01 and BF21 to BF10? Yes.1 
    # - Check style with lintR
    # - Think a better name for the condition object.
    # - What are the hypotheses that are going to be tested?
    # - Add plots.This could be a function.
    # - Add binary search algorithm.
    # - Run a simulation to know see the performance under various conditions.
    # - A>B vs B>A


# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed are not TRUE both.

# Test -------------------------------------------------------------------------
start.time <- Sys.time()
SSD_crt(eff.size = 0.5, n.datasets = 10, rho = 0.1, BF.thresh = 3, hypotheses = "vs.null",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE) #singularity
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.01, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
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

