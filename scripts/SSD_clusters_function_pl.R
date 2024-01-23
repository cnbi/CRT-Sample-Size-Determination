########## SSD NULL HYPOTHESIS AND SAVING EACH INCREASE #########################

# SSD_crt_null_plots <- function(eff_size, n1 = 15, n2 = 30, ndatasets = 1000, rho, BF_thresh,
#                                eta = 0.8, fixed = "n2", b_fract = 3, max = 1000) {
#   # Libraries ----
#   library(lme4)
#   #library(bain)
#   library(dplyr)
#   
#   # Warnings
#   if (is.numeric(c(eff_size, n1, n2, ndatasets, rho, BF_thresh, eta, b_fract)) == FALSE) stop("The arguments, wtih exception of fixed, must be numeric")
#   if (eff_size > 1) stop("Effect size must be standardized")
#   if (n2 %% 2 > 0) stop("Number of clusters must be even")
#   if (rho > 1) stop("The intraclass correlation must be standardized. Thus it cannot be more than 1")
#   if (eta > 1) stop("Probability of exceeding Bayes Factor threshold cannot be more than 1")
#   if (is.character(fixed) == FALSE) stop("Can only be a character indicating n1 or n2.")
#   if (fixed %in% c("n1", "n2") == FALSE) stop("Can only be a character indicating n1 or n2.")
#   if ((b_fract == round(b_fract)) == FALSE) stop("The fraction of information (b) must be integer")
#   
#   #Functions ----------------
#   source("data_generation.R")
#   source("small_functions.R")
#   source("print_results.R")
#   source("aafbf.R")
#   
#   # Starting values ----------------------------------------------------------
#   total_var <- 1
#   var_u0 <- rho * total_var       #Between-cluster variance
#   var_e <- total_var - var_u0     #Within-cluster variance
#   #iterations <- 1 DELETE THIS
#   eff_size0 <- 0                  #Effect size for null hypothesis
#   condition_met <- FALSE          #Indication we met the power criteria.
#   best_result <- FALSE            #Indication that we found the optimal value of sample size.
#   
#   
#   # Binary search start ------------------------------
#   if (fixed == "n1") {
#     low <- 6                   #lower bound
#   } else if (fixed == "n2") {
#     low <- 5                   #lower bound
#   }
#   high <- max                    #higher bound
#   
#   #Hypotheses -----------------------------------------
#   hypothesis1 <- "Intervention>Control"
#   null <- "Intervention=Control"
#   # hypoth <- paste(null, ";", hypothesis1)
#   
#   final_SSD <- vector(mode = "list", length = b_fract + 2)  # 480 = ((high-n2)/2)*b_fract
#   SSD_object <- vector(mode = "list", length = 321)   # CHANGE THIS
#   # Simulation and evaluation of condition  ----------------------------------
#   for (b in seq(b_fract)) {
#     high <- max
#     iteration <- 0
#     while (n2 < high + 1) {
#       # If H1 is true
#       data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
#                                             mean_interv = eff_size, b,
#                                             type = "equality"))
#       colnames(data_H1) <- c("BF.10", "BF.01",
#                              "PMP.0", "PMP.1")
#       # If H0 is true
#       data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
#                                             mean_interv = eff_size0, b,
#                                             type = "equality"))
#       colnames(data_H0) <- c("BF.10", "BF.01",
#                              "PMP.0", "PMP.1")
#       #Evaluation of condition ----
#       # Proportion
#       prop_BF01 <- length(which(data_H0[, "BF.01"] > BF_thresh)) / ndatasets
#       prop_BF10 <- length(which(data_H1[, "BF.10"] > BF_thresh)) / ndatasets
#       # save results
#       iteration <- iteration + 1
#       SSD_object[[iteration]] <- list("n1" = n1,
#                                       "n2" = n2,
#                                       "Proportion.BF01" = prop_BF01,
#                                       "Proportion.BF10" = prop_BF10,
#                                       "b.frac" = b,
#                                       "data_H0" = data_H0,
#                                       "data_H1" = data_H1)
#       # Increase
#       if (fixed == "n1") {
#         #Increase number of clusters
#         n2 = n2 + 2
#       } else if (fixed == "n2") {
#         n1 = n1 + 1
#       }
#       # Evaluation
#       # ifelse(prop_BF01 > eta & prop_BF10 > eta, condition_met <- TRUE, condition_met <- FALSE)
#       # actual <- min(prop_BF10, prop_BF01)
#       # Binary search algorithm -----
#       # if (condition_met == FALSE) {
#       #     print(c("Using cluster size:", n1, "and number of clusters:", n2,
#       #             "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10, "b:", b))
#       #     if (fixed == "n1") {      # We need to increase the sample size
#       #         low <- n2                         #lower bound
#       #         high <- high                      #higher bound
#       #         n2 <- round((low + high) / 2)     #point in the middle
#       #         ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
#       #         if (low + n2 == high * 2){          #when there is a roof effect
#       #             low <- n2
#       #             high <- max
#       #             n2 <- round((low + high) / 2)     #point in the middle
#       #         }
#       #     } else if (fixed == "n2") { # We need to increase the sample size
#       #         low <- n1                        #lower bound
#       #         high <- high                     #higher bound
#       #         n1 <- round((low + high) / 2)    #point in the middle
#       #         if (low + n1 == high * 2){         #when there is a roof effect
#       #             low <- n1
#       #             high <- previous_high
#       #             n1 <- round((low + high) / 2)    #point in the middle
#       #         }
#       # 
#       #     }
#       # } else if (condition_met == TRUE) {
#       #     print(c("Using cluster size:", n1,
#       #             "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
#       #             "low: ", low, "n2: ", n2, "high: ", high, "b:", b))
#       #     previous_high <- high
#       #     if (fixed == "n1") {
#       #         if (actual - eta < 0.1) { # Proportion is close enough
#       #             best_result == TRUE
#       #             break
#       #         } else if (previous == actual) { #If there is no change in the proportion and the lower
#       #             if (n2 - low == 2) {         #bound is close to the middle point.
#       #                 best_result == TRUE
#       #                 break
#       #             } else {   # Decreasing to find the optimal sample size
#       #                 low <- low
#       #                 high <- n2
#       #                 n2 <- round((low + high) / 2)   #point in the middle
#       #                 ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
#       #                 if (n2 < 30) warning("The number of groups is less than 30.
#       #                                      This could lead to problems in convergence and singularity.")
#       #                 print("Lowering") # Eliminate later
#       #             }
#       #         } else {       # Decreasing to find the optimal sample size
#       #             low <- low
#       #             high <- n2
#       #             n2 <- round((low + high) / 2)    #point in the middle
#       #             ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
#       #             if (n2 < 30) warning("The number of groups is less than 30.
#       #                                  This could lead to problems in convergence and singularity.")
#       #             print("Lowering") # Eliminate later
#       #         }
#       #     } else if (fixed == "n2") {
#       #         if (actual - eta < 0.1) { # Proportion close enough
#       #             best_result == TRUE
#       #             break
#       #         } else if (actual == previous) { # There is no change in the proportion and the lower bound
#       #             if (n1 - low == 1) {         #bound is close enough
#       #                 best_result == TRUE
#       #                 break
#       #             } else {                # Decreasing to find the optimal
#       #                 low <- low
#       #                 high <- n1
#       #                 n1 <- round((low + high) / 2)   #point in the middle
#       #                 print("Lowering") # Eliminate later
#       #             }
#       #         } else {                    # Decreasing to find the optimal
#       #             low <- low
#       #             high <- n1
#       #             n1 <- round((low + high) / 2)    #point in the middle
#       #             print("Lowering") # Eliminate later
#       #         }
#       #     }
#       # }
#       print(c("low:", low, "n2:", n2, "h:", high, "b:", b)) # Eliminate
#     }
#     final_SSD[[b]] <- SSD_object
#     print(c("Using cluster size: ", n1, "and number of clusters: ", n2,
#             "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10))
#   }
#   final_SSD[[b_fract + 1]] <- list(null, hypothesis1)
#   final_SSD[[b_fract + 2]] <- BF_thresh
#   
#   # Final output -----
#   print_results(final_SSD)
#   invisible(final_SSD)
# }
# 


#-----------------------------------------------------------------------------------------------
SSD_crt_plots <- function(eff_size, n1 = 15, n2 = 30, ndatasets = 1000, rho, BF_thresh,
                          eta = 0.8, fixed = "n2", b_fract = 3, max = 1000, batch_size = 100) {
    # Libraries ----
    library(lme4)
    library(dplyr)
    
    # Warnings
    if (is.numeric(c(eff_size, n1, n2, ndatasets, rho, BF_thresh, eta, b_fract, max, batch_size)) == FALSE) 
        stop("All arguments, except 'fixed', must be numeric")
    if (eff_size < 0) stop("The effect size must be a positive value ")
    if (n2 %% 2 > 0) stop("The number of clusters must be even")
    if (rho > 1) stop("The intraclass correlation must be standardized and cannot be larger than 1")
    if (rho < 0) stop("The intraclass correlation must be a positive value")
    if (eta > 1) stop("The probability of exceeding Bayes Factor threshold cannot be larger than 1")
    if (eta < 0) stop("The probability of exceeding Bayes Factor threshold must be a positive value")
    if (is.character(fixed) == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if (fixed %in% c("n1", "n2") == FALSE) stop("Fixed can only be a character indicating n1 or n2.")
    if ((b_fract == round(b_fract)) == FALSE) stop("The fraction of information (b) must be an integer")
    
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
    ultimate_sample_size <- FALSE            #Indication that we found the required sample size.
    
    
    # Binary search start ------------------------------------------------------
    if (fixed == "n1") {
        min_sample <- 6                     # Minimum cluster size
        low <- min_sample                   #lower bound
    } else if (fixed == "n2") {
        min_sample <- 5                     # Minimum number of clusters
        low <- min_sample                   #lower bound
    }
    high <- max                    #higher bound
    
    #Hypotheses ----------------------------------------------------------------
    hypothesis1 <- "Intervention>Control"
    null <- "Intervention=Control"
    final_SSD <- vector(mode = "list", length = b_fract)
    
    # Simulation of data and evaluation of condition  ----------------------------------
    low <- min_sample
    previous_eta <- 0
    previous_high <- 0
    high <- max
    # If H1 is true
    data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                          mean_interv = eff_size, b_fract, 
                                          type = "equality", batch_size = batch_size))
    colnames(data_H1) <- c("BF.10", "BF.01",
                           "PMP.0", "PMP.1")
    # If H0 is true
    data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                          mean_interv = eff_size0, b_fract,
                                          type = "equality", batch_size == batch_size))
    colnames(data_H0) <- c("BF.10", "BF.01",
                           "PMP.0", "PMP.1")
    print("Data generation check")
    #Evaluation of condition -------------------------------------------
    # Proportion
    prop_BF01 <- length(which(data_H0[, "BF.01"] > BF_thresh)) / ndatasets
    prop_BF10 <- length(which(data_H1[, "BF.10"] > BF_thresh)) / ndatasets
    # Save results
    SSD_object <- list("n1" = n1,
                       "n2" = n2,
                       "Proportion.BF01" = prop_BF01,
                       "Proportion.BF10" = prop_BF10,
                       "b.frac" = b,
                       "data_H0" = data_H0,
                       "data_H1" = data_H1)
    final_SSD <- SSD_object
    rm(data_H0, data_H1)
    
    # final_SSD[[b_fract + 1]] <- list(null, hypothesis1)
    # final_SSD[[b_fract + 2]] <- BF_thresh
    
    # Final output -----
    #print_results(final_SSD)
    invisible(final_SSD)
    #     # Binary search algorithm ------------------------------------------
    #     if (condition_met == FALSE) {
    #         print(c("Using cluster size:", n1, "and number of clusters:", n2,
    #                 "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10, "b:", b))
    #         if (fixed == "n1") {
    #             # Increase the number of clusters since eta is too small
    #             low <- n2                         #lower bound
    #             high <- high                      #higher bound
    #             n2 <- round((low + high) / 2)     #point in the middle
    #             ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1) # Ensure number of clusters is even
    #             
    #             # Adjust higher bound when there is a ceiling effect
    #             if (low + n2 == high * 2) {
    #                 low <- n2                         #lower bound
    #                 high <- max                       #higher bound
    #                 n2 <- round((low + high) / 2)     #point in the middle
    #             }
    #         } else if (fixed == "n2") {
    #             # Increase the cluster sizes since eta is too small
    #             low <- n1                        #lower bound
    #             high <- high                     #higher bound
    #             n1 <- round((low + high) / 2)    #point in the middle
    #             
    #             # Adjust higher bound when there is a ceiling effect
    #             if (low + n1 == high * 2) {
    #                 low <- n1                        #lower bound
    #                 #Set the higher bound based on the previous high or the maximum
    #                 if (previous_high > 0) {
    #                     high <- previous_high
    #                 } else {
    #                     high <- max
    #                 }
    #                 n1 <- round((low + high) / 2)    #point in the middle
    #             }
    #         }
    #     } else if (condition_met == TRUE) {
    #         print(c("Using cluster size:", n1,
    #                 "and number of clusters:", n2,
    #                 "prop_BF01: ", prop_BF01, "prop_BF10: ", prop_BF10,
    #                 "low: ", low, "high: ", high, "b:", b))
    #         previous_high <- high
    #         
    #         if (fixed == "n1") {
    #             # Eta is close enough to the desired eta
    #             if (actual_eta - eta < 0.1) {
    #                 if (previous_eta > eta) {
    #                     # Decreasing to find the ultimate number of clusters
    #                     low <- low                         #lower bound
    #                     high <- n2                         #higher bound
    #                     n2 <- round((low + high) / 2)      #point in the middle
    #                     ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
    #                     if (n2 < 30) warning("The number of groups is less than 30.
    #                                          This may cause problems in convergence and singularity.")
    #                     print("Lowering") # Eliminate later
    #                 }else {
    #                     ultimate_sample_size == TRUE
    #                     break
    #                 }
    #             } else if (previous_eta == actual_eta) {
    #                 # If there is no change in eta and the lower bound is close to the middle point
    #                 if (n2 - low == 2) {
    #                     ultimate_sample_size == TRUE
    #                     break
    #                     
    #                 } else {   
    #                     # Decreasing to find the ultimate number of clusters
    #                     low <- low                         #lower bound
    #                     high <- n2                         #higher bound
    #                     n2 <- round((low + high) / 2)      #point in the middle
    #                     ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
    #                     if (n2 < 30) warning("The number of groups is less than 30.
    #                                          This may cause problems in convergence and singularity.")
    #                     print("Lowering") # Eliminate later
    #                 }
    #             } else {
    #                 # Decreasing to find the ultimate number of clusters
    #                 low <- low                         #lower bound
    #                 high <- n2                         #higher bound
    #                 n2 <- round((low + high) / 2)      #point in the middle
    #                 ifelse(n2 %% 2 == 0, n2 <- n2, n2 <- n2 + 1)
    #                 if (n2 < 30) warning("The number of groups is less than 30.
    #                                      This may cause problems in convergence and singularity.")
    #                 print("Lowering") # Eliminate later
    #             }
    #         } else if (fixed == "n2") {
    #             # Eta is close enough to the desired eta
    #             if (actual_eta - eta < 0.1) {
    #                 ultimate_sample_size == TRUE
    #                 break
    #                 
    #             } else if (actual_eta == previous_eta) {
    #                 # If there is no change in eta and the lower bound is close to the middle point
    #                 if (n1 - low == 1) {
    #                     ultimate_sample_size == TRUE
    #                     break
    #                     
    #                 } else {
    #                     # Decreasing the cluster size to find the ultimate sample size
    #                     low <- low                         #lower bound
    #                     high <- n1                         #higher bound
    #                     n1 <- round((low + high) / 2)      #point in the middle
    #                     print("Lowering") # Eliminate later
    #                 }
    #             } else {
    #                 # Decreasing the cluster size to find the ultimate sample size
    #                 low <- low                         #lower bound
    #                 high <- n1                         #higher bound
    #                 n1 <- round((low + high) / 2)      #point in the middle
    #                 print("Lowering") # Eliminate later
    #             }
    #         }
    #     }
    #     previous_eta <- min(prop_BF10, prop_BF01)
    #     print(c("low:", low, "n2:", n2, "n1:", n1, "h:", high, "b:", b)) # Eliminate
    #     
    #     # If the sample size reaches the maximum
    #     if (n2 == max) {
    #         break
    #     } else if (n1 == max) {
    #         break
    #     }
    # }

}
