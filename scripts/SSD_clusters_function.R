########### Sample Size Determination for Cluster Randomized Trials ############

#"@title Sample Size Determination for Cluster Randomized Trials
#"@description
#"@arguments
## eff_size: Numeric. Effect size
## n1: Numeric. Cluster size
## n2: Numeric. Number of clusters in each treatment condition.
## ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size.
## rho: Numeric. Intraclass correlation
## BF_thresh: Numeric. Value of the Bayes factor that is going to be the threshold.
## eta: Numeric. Probability of exceeding Bayes Factor threshold.
## fixed: Character. Indicating which sample is fixed (n1 or n2)
## b_fract: Numeric. Fraction of information used to specify the prior distribution.


SSD_crt_null <- function(eff_size, n1 = 15, n2 = 30, ndatasets = 1000, rho, BF_thresh,
                         eta = 0.8, fixed = "n2", b_fract = 3) {
    # Libraries ----
    library(lme4)
    #library(bain)
    library(dplyr)

    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, rho, BF_thresh, eta, b_fract)) == FALSE) stop("The arguments, wtih exception of fixed, must be numeric")
    if (eff_size > 1) stop("Effect size must be standardized")
    if (n2 %% 2 > 0) stop("Number of clusters must be even")
    if (rho > 1) stop("The intraclass correlation must be standardized. Thus it cannot be more than 1")
    if (eta > 1) stop("Probability of exceeding Bayes Factor threshold cannot be more than 1")
    if (is.character(fixed) == FALSE) stop("Can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Can only be a character indicating n1 or n2.")
    if ((b_fract == round(b_fract)) == FALSE) stop("The fraction of information (b) must be integer")

    #Functions ----------------
    source("data_generation.R")
    source("small_functions.R")
    source("print_results.R")
    source("aafbf.R")

    # Starting values ----------------------------------------------------------
    total_var <- 1
    var_u0 <- rho * total_var       #Between-cluster variance
    var_e <- total_var - var_u0     #Within-cluster variance
    #iterations <- 1 DELETE THIS
    eff_size0 <- 0                  #Effect size for null hypothesis
    condition_met <- FALSE          #Indication we met the power criteria.
    best_result <- FALSE            #Indication that we found the optimal value of sample size.


    # Binary search start ------------------------------
    if (fixed == "n1") {
        low <- 6                   #lower bound
    } else if (fixed == "n2") {
        low <- 5                   #lower bound
    }
    high <- 1000                    #higher bound

    #Hypotheses -----------------------------------------
    # hypothesis1 <- "Dintervention>Dcontrol"
    # null <- "Dintervention=Dcontrol"
    # hypoth <- paste(null, ";", hypothesis1)
    
    final_SSD <- vector(mode = "list", length = b_fract)

    # Simulation and evaluation of condition  ----------------------------------
    for (b in seq(b_fract)) {
        previous <- 0
        while (best_result == FALSE) {
            # If H1 is true
            data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                                  mean_interv = eff_size, b, 
                                                  type = "equality"))
            colnames(data_H1) <- c("BF.10", "BF.01",
                                   "PMP.0", "PMP.1")
            # If H0 is true
            data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                                  mean_interv = eff_size0, b,
                                                  type = "equality"))
            colnames(data_H0) <- c("BF.10", "BF.01",
                                   "PMP.0", "PMP.1")
            #Evaluation of condition ----
            # Proportion
            prop_BF01 <- length(which(data_H0[, "BF.01"] > BF_thresh)) / ndatasets
            prop_BF10 <- length(which(data_H1[, "BF.10"] > BF_thresh)) / ndatasets
            # Evaluation
            ifelse(prop_BF01 > eta & prop_BF10 > eta, condition_met <- TRUE, condition_met <- FALSE)
            actual <- min(prop_BF10, prop_BF01)
            # Binary search algorithm -----
            if (condition_met == FALSE) {
                print(c("Using cluster size:", n1, "and number of clusters:", n2,
                        "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10, "b:", b))
                if (fixed == "n1") {      # We need to increase the sample size
                    low <- n2                         #lower bound
                    high <- high                      #higher bound
                    n2 <- round((low + high) / 2)     #point in the middle
                    ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                } else if (fixed == "n2") { # We need to increase the sample size
                    low <- n1                        #lower bound
                    high <- high                     #higher bound
                    n1 <- round((low + high) / 2)    #point in the middle
                }
            } else if (condition_met == TRUE) {
                print(c("Using cluster size:", n1,
                        "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
                        "low: ", low, "n2: ", n2, "high: ", high, "b:", b))
                if (fixed == "n1") {
                    if (actual - eta < 0.1) { # Proportion is close enough
                        best_result == TRUE
                        break
                    } else if (previous == actual) { #If there is no change in the proportion and the lower
                        if (n2 - low == 2) {         #bound is close to the middle point.
                            best_result == TRUE
                            break
                        } else {   # Decreasing to find the optimal sample size
                            low <- low
                            high <- n2
                            n2 <- round((low + high) / 2)   #point in the middle
                            ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                            if (n2 < 30) warning("The number of groups is less than 30.
                                                 This could lead to problems in convergence and singularity.")
                            print("Lowering") # Eliminate later
                        }
                    } else {       # Decreasing to find the optimal sample size
                        low <- low
                        high <- n2
                        n2 <- round((low + high) / 2)    #point in the middle
                        ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
                        if (n2 < 30) warning("The number of groups is less than 30.
                                             This could lead to problems in convergence and singularity.")
                        print("Lowering") # Eliminate later
                    }
                } else if (fixed == "n2") {
                    if (actual - eta < 0.1) { # Proportion close enough
                        best_result == TRUE
                        break
                    } else if (actual == previous) { # There is no change in the proportion and the lower bound
                        if (n1 - low == 1) {         #bound is close enough
                            best_result == TRUE
                            break
                        } else {                # Decreasing to find the optimal
                            low <- low
                            high <- n1
                            n1 <- round((low + high) / 2)   #point in the middle
                            print("Lowering") # Eliminate later
                        }
                    } else {                    # Decreasing to find the optimal
                        low <- low
                        high <- n1
                        n1 <- round((low + high) / 2)    #point in the middle
                        print("Lowering") # Eliminate later
                    }
                }
            }
            previous <- min(prop_BF10, prop_BF01)
            print(c("low:", low, "n2:", n2, "h:", high, "b:", b)) # Eliminate
            if (n2 == 1000) {
                break
            } else if (n1 == 1000) {
                break
            }
        }
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF01" = prop_BF01,
                           "Proportion.BF10" = prop_BF10,
                           "b.frac" = b,
                           "data_H0" = data_H0,
                           "data_H1" = data_H1)
        final_SSD[[b]] <- SSD_object
    }
    final_SSD[[b_fract + 1]] <- list(null, hypothesis1)
    final_SSD[[b_fract + 2]] <- BF_thresh

    # Final output -----
    print_results(final_SSD)
    invisible(final_SSD)
}



#TODO:



# Test -------------------------------------------------------------------------
start.time <- Sys.time()
nulla <- SSD_crt_null(eff_size = 0.5, ndatasets = 10, rho = 0.1, BF_thresh = 3, fixed = "n1",
                      b_fract = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff_size = 0.4, ndatasets = 15, rho = 0.05, BF_thresh = 3, fixed = "n1", b_fract = 2) #singularity
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff_size = 0.4, ndatasets = 15, rho = 0.01, BF_thresh = 3, fixed = "n1", b_fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # # This is weird, it seems like this needs less number of clusters.
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff_size = 0.4, ndatasets = 15, rho = 0.1, BF_thresh = 4, fixed = "n1", b_fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
#
# start.time <- Sys.time()
# try <- SSD_crt_null(eff_size = 0.2, ndatasets = 20, rho = 0.1, BF_thresh = 5,fixed = "n1",
#                     b_fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # 
# start.time <- Sys.time()
# SSD_crt_null(eff_size = 0.2, ndatasets = 20, rho = 0.05, BF_thresh = 3, fixed = "n2", b_fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# Reach maximum value
# 
# start.time <- Sys.time()
# SSD_crt_null(eff_size = 0.8, ndatasets = 20, rho = 0.1, BF_thresh = 3, fixed = "n2", b_fract = 3)
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken
# # Again maximum value
