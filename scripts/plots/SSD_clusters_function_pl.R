########## SSD NULL HYPOTHESIS AND SAVING EACH INCREASE #########################

#"@arguments
## eff_size: Numeric. Effect size
## n1: Numeric. Cluster size
## n2: Numeric. Number of clusters per treatment condition.
## ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size.
## rho: Numeric. Intraclass correlation
## BF_thresh: Numeric. Value of the Bayes factor that is going to be the threshold.
## eta: Numeric. Probability of exceeding Bayes Factor threshold.
## fixed: Character. Indicating which sample is fixed (n1 or n2)
## b_fract: Numeric. Fraction of information used to specify the prior distribution.
## max: Maximum sample size.
## bacth_size: This parameter determines the size of batches used during the fitting of the multilevel model.

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
    
    # Functions ----------------
    source("data_generation.R")
    source("small_functions.R")
    source("print_results.R")
    source("aafbf.R")
    
    # Starting values ----------------------------------------------------------
    total_var <- 1
    var_u0 <- rho * total_var       #Between-cluster variance
    var_e <- total_var - var_u0     #Within-cluster variance
    eff_size0 <- 0                  #Effect size for null hypothesis
    
    # Results ----------------------------------------------------------------
    final_SSD <- vector(mode = "list", length = b_fract)
    
    # Simulation of data and evaluation of condition  --------------------------
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
                       "b.frac" = b_fract,
                       "data_H0" = data_H0,
                       "data_H1" = data_H1)
    final_SSD <- SSD_object
    rm(data_H0, data_H1)
    
    invisible(final_SSD)

}
