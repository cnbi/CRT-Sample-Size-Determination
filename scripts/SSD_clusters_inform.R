############# Sample Size Determination for Cluster Randomized Trials evaluating  #############
############################### informative hypotheses ######################################

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description 
#'@arguments
## eff.size: Numeric. Effect size
## n1: Numeric. Cluster size
## n2: Numeric. Number of clusters in each treatment condition.
## n.datasets: Numeric. Number of data sets that the user wants to generate to determine the sample size.
## rho: Numeric. Intraclass correlation
## BF.thresh: Numeric. Value of the Bayes factor that is going to be the threshold.
## eta: Numeric. Probability of exceeding Bayes Factor threshold.
## fixed: Character. Indicating which sample is fixed (n1 or n2)

SSD_crt_inform <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh,
                           eta = 0.8, fixed = 'n2') {
    # Libraries -----------------
    library(lme4)
    #library(bain)
    library(dplyr)
    
    # Warnings
    if (is.numeric(c(eff.size, n1, n2, n.datasets, rho, BF.thresh, eta)) == FALSE) stop("The arguments, wtih exception of fixed, must be numeric") 
    if (eff.size > 1) stop("Effect size must be standardized")
    if (n2 %% 2 > 0) stop("Number of clusters must be even")
    if (rho > 1) stop("The intraclass correlation must be standardized. Thus it cannot be more than 1")
    if (eta > 1) stop("Probability of exceeding Bayes Factor threshold cannot be more than 1")
    if (is.character(fixed) == FALSE) stop("Can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Can only be a character indicating n1 or n2.")
    
    # Functions
    source('data_generation.R')
    source("small_functions.R")
    source("print_results.R")
    source("aafbf.R")
    
    # Starting values ------------
    total.var <- 1
    var.u0 <- rho * total.var       #Between-cluster variance
    var.e <- total.var - var.u0     #Within-cluster variance
    iterations <- 1
    condition_met <- FALSE          #Indicator of the fulfillment of the power criteria.
    best_result <- FALSE            #Indicator if we find the optimal value of sample size.
    
    # Binary search start
    if (fixed == 'n1') {
        low <- 6                   #lower bound
    } else if (fixed == 'n2') {
        low <- 5                   #lower bound
    }
    high <- 1000                    #higher bound
    
    #Hypotheses
    hypothesis1 <- "Dintervention>Dcontrol"
    hypothesis2 <- "Dintervention<Dcontrol"
    #hypoth <- paste(hypothesis1, ";", hypothesis2)
    
    previous <- 0
    
    # Simulation ---------------
    while (best_result == FALSE) {
        # If H1 is true
        data_crt <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, 
                                               eff.size, b = 1, type = "inequalities"))
        colnames(data_crt) <- c('BF.12', 'BF.21', 'PMP.1', 'PMP.2')
        # Evaluation of condition----
        # Proportion
        # browser()
        prop.BF12 <- length(which(data_crt[, 'BF.12'] > BF.thresh)) / n.datasets 
        # Evaluation
        ifelse(prop.BF12 > eta, condition_met <- TRUE, condition_met <- FALSE) #Maybe I have to change this part
        # browser()
        # Binary search algorithm----
        if (condition_met == FALSE) {
            print(c("Using cluster size:", n1, "and number of clusters:", n2,
                  "prop.BF12: ", prop.BF12))
            if (fixed == 'n1') {
                low <- n2                        #lower bound
                high <- high                     #higher bound
                n2 <- round((low + high) / 2)    #point in the middle
                ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
            } else if (fixed == 'n2') {
                low <- n1                        #lower bound
                high <- high                     #higher bound
                n1 <- round((low + high) / 2)    #point in the middle
            }
        } else if (condition_met == TRUE) {
            print(c("Using cluster size:", n1, "and number of clusters:", n2,
                  "prop.BF12:", prop.BF12))
            if (fixed == 'n1') {
                if (prop.BF12 - eta < 0.1) { # Proportion is close enough
                    best_result == TRUE
                    break
                } else if (previous == prop.BF12) {
                    if (n2 - low == 2) {
                        best_result == TRUE
                        print("Found it!") # Eliminate later
                        break
                    } else { # Decreasing
                        low <- low                         #lower bound
                        high <- n2                       #higher bound
                        n2 <- round((low + high) / 2)    #point in the middle
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                        print("Lowering") # Eliminate later
                        if (n2 < 30) warning("The number of groups is less than 30. This could lead to problems in convergence and singularity.")
                    }
                } else {   #Decreasing to find the optimal
                    low <- low                       #lower bound
                    high <- n2                       #higher bound
                    n2 <- round((low + high) / 2)    #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                    if (n2 < 30) warning("The number of groups is less than 30. This could lead to problems in convergence and singularity.")
                    print("Lowering") # Eliminate later
                }
            } else if (fixed == 'n2') {
                if (prop.BF12 - eta < 0.1) { # Proportion close enough
                    best_result == TRUE
                    break
                } else if (prop.BF12 == previous) { # There is no change in the proportion and the lower bound
                    if (n1 - low == 1) {
                        best_result == TRUE
                        break
                    } else {
                        low <- low                         #lower bound
                        high <- n1                       #higher bound
                        n1 <- round((low + high) / 2)    #point in the middle
                        print("Lowering") # Eliminate later
                    }
                } else {
                    low <- low                       #lower bound
                    high <- n1                       #higher bound
                    n1 <- round((low + high) / 2)    #point in the middle
                    print("Lowering") # Eliminate later
                }
            }
        }
        previous <- prop.BF12
        # Stop because the number of clusters is crazy, is not plausible
        if (n2 == 1000) {
            break
        } else if (n1 == 1000) {
            break
        }
        
    }
    # Final output -----------
    SSD_object <- list("n1" = n1,
                       "n2" = n2,
                       "Eta" = prop.BF12,
                       "data" = data_crt,
                       "hypotheses" = list(hypothesis1, hypothesis2),
                       "BF.threshold" = BF.thresh,
                       "evaluation" = "Inequalities")
    
    print_results(SSD_object)
    invisible(SSD_object)
}


#TODO:
# - Check style with lintR
# - Add plots.This could be a function.
# - Run a simulation to know see the performance under various conditions.


# Test -------------------------------------------------------------------------
# start.time <- Sys.time()
# a <- SSD_crt_inform(eff.size = 0.5, n.datasets = 100, rho = 0.1, BF.thresh = 3, fixed = 'n1')
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# start.time <- Sys.time()
# SSD_crt_inform(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, fixed = "n1") #singularity
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# 
# start.time <- Sys.time()
# SSD_crt_inform(eff.size = 0.4, n.datasets = 15, rho = 0.01, BF.thresh = 6, fixed = "n1")
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # This is weird, it seems like this needs less number of clusters.
# 
# start.time <- Sys.time()
# SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
#         n1.fixed = TRUE, n2.fixed = FALSE)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# start.time <- Sys.time()
# SSD_crt(eff.size = 0.2, n.datasets = 20, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
#         n1.fixed = TRUE, n2.fixed = FALSE)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# start.time <- Sys.time()
# SSD_crt(eff.size = 0.2, n.datasets = 20, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
#         n1.fixed = TRUE, n2.fixed = FALSE)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
# start.time <- Sys.time()
# SSD_crt(eff.size = 0.8, n.datasets = 20, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
#         n1.fixed = TRUE, n2.fixed = FALSE)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# 
