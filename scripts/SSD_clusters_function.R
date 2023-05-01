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


SSD_crt <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh, 
                    eta = 0.8, hypothesis, n1.fixed = TRUE,
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

    #Hypothesis
    if (hypothesis == "interv.bigger") {
        hypothesis <- "Dintervention>Dcontrol"
    } else if (hypothesis == "interv.smaller") {
        hypothesis <- "Dintervention<Dcontrol" 
    }
    null <- "Dintervention=Dcontrol"
    hypoth <- paste(null, ";", hypothesis)
    # Starting values
    total.var <- 1
    mean.ctrl <- 0
    var.u0 <- rho * total.var
    var.e <- total.var - var.u0
    iterations <- 1
    eff.size0 <- 0
    
    # Names of table with results
    names.table <- c('Dcontrol', 'Dintervention', 'BF.12', 'BF.21', 'BFc.1', 'BFc.2', 'PMP.1',
                     'PMP.2', 'PMPc.1', 'PMPc.2', 'PMPc.c')
    condition <- FALSE #condition fulfillment indicator
    
    while (condition == FALSE) {
        # If H1 is true
        data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho, 
                                              eff.size, hypoth, mean.ctrl))
        print("Yay data ready!")
        results.H1 <- do.call(extract_results, list(n.datasets, data.H1, names.table = names.table))
        print("Yay extract results is ready!")
        # If H0 is true
        data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, rho,
                                              eff.size = eff.size0, hypoth, mean.ctrl))
        results.H0 <-  do.call(extract_results, list(n.datasets, data.H0, names.table = names.table))
        #Evaluation of condition
        condition <- eval_thresh(results.H0 = results.H0, results.H1 = results.H1, 
                                 BF.thresh = BF.thresh, n.datasets = n.datasets, condition = condition, eta = eta)
        print("Yay evaluation works!")
        # If condition =FALSE then increase the sample. Which n increase?
        if (condition == FALSE){
            if (n1.fixed == TRUE) {
                n2 = n2 + 2
            } else if (n2.fixed == TRUE) {
                n1 = n1 + 1
            }
        }
        iterations <- iterations + 1
        print(c("n1: ", n1, " n2: ", n2))
        if (n2 == 1000) {
            break
        }
    }
    
    # output
    print(c("Cluster size: ", n1, "Number of clusters: ", n2))
    
}


#TODO:
    # - Add ICC calculation
    # - Think: If n1.fixed and n2.fixed are FALSE, then what?
    # - Think: How is going to be the output?
    # - Think: Should I change BF12 to BF01 and BF21 to BF10?
    # - Check style with lintR
    # - Test evaluation function.
    # - Think a better name for the condition object.
    # - Hypothesis should be entered by researcher. Change this in data_generation function.
    # - What are the hypotheses that are going to be tested?
    # - Add plots.This could be a function.
    # - Add binary search algorithm.
    # - Right now only woks for the hypothesis Dintervention > Dcontrol. Generalize a little more.
    # - Run a simulation to know see the performance under various conditions.



# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed are not TRUE both.

# Test -------------------------------------------------------------------------
start.time <- Sys.time()
SSD_crt(eff.size = 0.5, n.datasets = 10, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
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

