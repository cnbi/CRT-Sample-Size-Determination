############# Sample Size Determination for Cluster Randomized Trials #############

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
## b.fract: Numeric. Fraction of information used to specify the prior distribution.


SSD_crt_null <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh,
                         eta = 0.8, fixed = 'n2', b.fract = 3) {
    # Libraries ----
    library(lme4)
    library(bain)
    library(dplyr)
    
    # Warnings
    if (is.numeric(c(eff.size, n1, n2, n.datasets, rho, BF.thresh, eta, b.fract)) == FALSE) stop("The arguments, wtih exception of fixed, must be numeric") 
    #if (is.numeric(eff.size) == FALSE) stop("The effect size must be numeric")
    if (eff.size > 1) stop("Effect size must be standardized")
    #if (is.numeric(n1) == FALSE) stop("Cluster size must be numeric")
    #if (is.numeric(n2) == FALSE) stop("Number of clusters must be numeric")
    if (n2 %% 2 > 0) stop("Number of clusters must be even")
    #if (is.numeric(n.datasets) == FALSE) stop("The number of data sets must be numeric")
    if (rho > 1) stop("The intraclass correlation must be standardized. Thus it cannot be more than 1")
    #if (is.numeric(rho) == FALSE) stop("The intraclass correlation must be numeric")
    #if (is.numeric(BF.thresh) == FALSE) stop("The threshold must be numeric")
    #if (is.numeric(eta) == FALSE) stop("The probability of exceeding Bayes Factor must be numeric")
    if (eta > 1) stop("Probability of exceeding Bayes Factor threshold cannot be more than 1")
    if (is.character(fixed) == FALSE) stop("Can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Can only be a character indicating n1 or n2.")
    #if (is.numeric(b.fract) == FALSE) stop("The fraction b must be numeric")
    if ((b.fract == round(b.fract)) == FALSE) stop("The fraction of information (b) must be integer")
    
    #Functions ----
    source('data_generation.R')
    source("small_functions.R")
    source("print_results.R")
    
    # Starting values -----
    total.var <- 1
    var.u0 <- rho * total.var       #Between-cluster variance
    var.e <- total.var - var.u0     #Within-cluster variance
    iterations <- 1
    eff.size0 <- 0                  #Effect size for null hypothesis
    condition_met <- FALSE          #Indicator of the fulfillment of the power criteria.
    best_result <- FALSE            #Indicator if we find the optimal value of sample size.
    
    
    # Binary search start ----
    if (fixed == 'n1') {
        low <- 6                   #lower bound
    } else if (fixed == 'n2') {
        low <- 5                   #lower bound
    }
    high <- 1000                    #higher bound
    
    #Hypotheses -----
    hypothesis1 <- "Dintervention>Dcontrol"
    null <- "Dintervention=Dcontrol"
    hypoth <- paste(null, ";", hypothesis1)
    final_SSD <- vector(mode = "list", length = b.fract)
    
    # Simulation and evaluation of condition  -----
    for (b in seq(b.fract)) {
        previous <- 0
        previous_n2 <- 0
        previous_high <- 0
        previous_n1 <- 0
        while (best_result == FALSE) {
            # If H1 is true
            data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e,
                                                  mean.interv = eff.size, hypoth, b))
            colnames(data.H1) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1', 'marker')
            marker_H1 <- sum(data.H1[, 'marker'])       #How many matrices were singular
            # If H0 is true
            data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e,
                                                  mean.interv = eff.size0, hypoth, b))
            colnames(data.H0) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1', 'marker')
            marker_H0 <- sum(data.H0[, 'marker'])
            
            #Evaluation of condition ----
            # Proportion
            prop.BF01 <- length(which(data.H0[, 'BF.01'] > BF.thresh)) / n.datasets
            prop.BF10 <- length(which(data.H1[, 'BF.10'] > BF.thresh)) / n.datasets
            # Evaluation
            ifelse(prop.BF01 > eta & prop.BF10 > eta, condition_met <- TRUE, condition_met <- FALSE)
            actual <- min(prop.BF10, prop.BF01)
            
            # Binary search algorithm -----
            if (condition_met == FALSE) {
                print(c("Using cluster size:", n1, "and number of clusters:", n2,
                        "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10))
                if (fixed == 'n1') {                        # We need to increase
                    # if (actual == previous) {
                    #     low <- n2
                    #     high <- previous_n2
                    #     n2 <- round(( low + high) / 2)    #point in the middle
                    #     ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                    # } else {
                    #     low <- n2                         #lower bound
                    #     high <- high                      #higher bound
                    #     n2 <- round(( low + high) / 2)    #point in the middle
                    #     ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                    # }
                    low <- n2                         #lower bound
                    high <- high                      #higher bound
                    n2 <- round(( low + high) / 2)    #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                } else if (fixed == 'n2') {
                    # if (actual == previous) {
                    #     low <- n1
                    #     high <- previous_n1
                    #     n1 <- round((low + high) / 2)
                    # } else {
                    #     low <- n1                        #lower bound
                    #     high <- high                     #higher bound
                    #     n1 <- round((low + high) / 2)    #point in the middle
                    # }
                    low <- n1                        #lower bound
                    high <- high                     #higher bound
                    n1 <- round((low + high) / 2)    #point in the middle
                }
            } else if (condition_met == TRUE) {
                print(c("Using cluster size:", n1, 
                        "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10, 
                        "low: ", low, "n2: ", n2, "high: ", high))
                previous_high <- high
                previous_low <- low
                if (fixed == 'n1') {
                    if (actual - eta < 0.1) { # If the initial sample size produce a proportion close enough
                        best_result == TRUE
                        break
                    } else if (previous == actual) {
                        if (n2 - low == 2) {
                            best_result == TRUE
                            break
                        } else {
                            low <- low
                            high <- n2
                            n2 <- round((low + high) / 2)    #point in the middle
                            ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                            previous_n2 <- n2
                            if (n2 < 30) warning("The number of groups is less than 30. This could lead to problems in convergence and singularity.")
                            print("Lowering") # Eliminate later
                        }
                    } else {                  # We need to decrease to find the optimal sample size
                        low <- low
                        high <- n2
                        n2 <- round((low + high) / 2)    #point in the middle
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                        previous_n2 <- n2
                        if (n2 < 30) warning("The number of groups is less than 30. This could lead to problems in convergence and singularity.")
                        print("Lowering") # Eliminate later
                    }
                } else if (fixed == 'n2') {
                    if (actual - eta < 0.1) { # If the initial sample size produce a proportion close enough
                        best_result == TRUE
                        break
                    } else if (actual == previous) { # There is no change because we found the optimal
                        if (n1 - low == 1) {
                            best_result == TRUE
                            break
                        } else {
                            low <- low
                            high <- n1
                            n1 <- round((low + high) / 2)    #point in the middle
                            previous_n1 <- n1
                            print("Lowering") # Eliminate later
                        }
                    } else {
                        low <- low
                        high <- n1
                        n1 <- round((low + high) / 2)    #point in the middle
                        previous_n1 <- n1
                        print("Lowering") # Eliminate later
                    }
                }
            }
            previous <- min(prop.BF10, prop.BF01)
            print(c("low:", low, "n2:", n2, "h:", high, "b:", b)) # Eliminate later
            if (n2 == 1000) {
                break
            } else if (n1 == 1000) {
                break 
            }
        }
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF01" = prop.BF01,
                           "Proportion.BF10" = prop.BF10,
                           "b.frac" = b,
                           "data.H0" = data.H0,
                           "data.H1" = data.H1,
                           "marker.H0" = marker_H0,
                           "marker.H1" = marker_H1)
        final_SSD[[b]] <- SSD_object
    }
    final_SSD[[b.fract + 1]] <- list(null, hypothesis1)
    final_SSD[[b.fract + 2]] <- BF.thresh
    
    # Final output -----
    print_results(final_SSD)
    invisible(final_SSD)
}




#TODO:
# - Check style with lintR


# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed are not TRUE both.

# Test -------------------------------------------------------------------------
# start.time <- Sys.time()
# nulla <- SSD_crt_null(eff.size = 0.5, n.datasets = 10, rho = 0.1, BF.thresh = 3, fixed = "n1",
#                       b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, fixed = "n1", b.fract = 2) #singularity
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff.size = 0.4, n.datasets = 15, rho = 0.01, BF.thresh = 3, fixed = "n1", b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # # This is weird, it seems like this needs less number of clusters.
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff.size = 0.4, n.datasets = 15, rho = 0.1, BF.thresh = 4, fixed = "n1", b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
#
# start.time <- Sys.time()
# try <- SSD_crt_null(eff.size = 0.2, n.datasets = 20, rho = 0.1, BF.thresh = 5,fixed = "n1",
#                     b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff.size = 0.2, n.datasets = 20, rho = 0.05, BF.thresh = 3, fixed = "n2", b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # Reach maximum value
# 
# start.time <- Sys.time()
# SSD_crt_null(eff.size = 0.8, n.datasets = 20, rho = 0.1, BF.thresh = 3, fixed = "n2", b.fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # Again maximum value
