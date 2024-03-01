########################### Calculate Bayesian Power with Equality Constraints #########################

#"@title Sample Size Determination for Cluster Randomized Trials
#"@description
#"@arguments
## eff_size: Numeric. Effect size
## n1: Numeric. Cluster size
## n2: Numeric. Number of clusters per treatment condition
## ndatasets: Numeric. Number of data sets that the user wants to generate to determine the sample size
## rho: Numeric. Intraclass correlation
## BF_thresh: Numeric. Value of the Bayes factor that is going to be the threshold
## eta: Numeric. Probability of exceeding Bayes Factor threshold
## b_fract: Numeric. Fraction of information used to specify the prior distribution
## bacth_size: This parameter determines the size of batches used during the fitting of the multilevel model.


get_power_CRT <- function(eff_size, n1, n2, ndatasets = 1000, rho, BF_thresh,
                          eta = 0.8, b_fract = 3, batch_size = 100) {
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
    eff_size0 <- 0                  #Effect size for null hypothesis
    results_H0 <- matrix(NA, nrow = ndatasets, ncol = 4)
    results_H1 <- matrix(NA, nrow = ndatasets, ncol = 4)
    
    #Hypotheses ----------------------------------------------------------------
    hypothesis1 <- "Intervention>Control"
    null <- "Intervention=Control"
    final_SSD <- vector(mode = "list", length = b_fract)
    type <- "Equality"
    b <- 1
    singular_warn <- 0
    
    # Simulation of data and evaluation of condition  ----------------------------------
    # If H1 is true
    data_H1 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                          mean_interv = eff_size,
                                          batch_size = batch_size))
    
    # If H0 is true
    data_H0 <- do.call(gen_CRT_data, list(ndatasets, n1, n2, var_u0, var_e,
                                          mean_interv = eff_size0,
                                          batch_size = batch_size))
    
    while (b < (b_fract + 1)) {
        #Approximated adjusted fractional Bayes factors------------------------------
        n_eff_H1 <- ((n1 * n2) / (1 + (n1 - 1) * data_H1$rho_data)) / 2
        output_AAFBF_H1 <- Map(calc_aafbf, type, data_H1$estimates, data_H1$cov_list, list(b), n_eff_H1)
        
        n_eff_H0 <- ((n1 * n2) / (1 + (n1 - 1) * data_H0$rho_data)) / 2
        output_AAFBF_H0 <- Map(calc_aafbf, type, data_H0$estimates, data_H0$cov_list, list(b), n_eff_H0)
        
        # Results ---------------------------------------------------------------------
        results_H1[, 1] <- unlist(lapply(output_AAFBF_H1, extract_res, 1)) # Bayes factor H1vsH0
        results_H1[, 2] <- unlist(lapply(output_AAFBF_H1, extract_res, 4)) #posterior model probabilities of H1
        results_H1[, 3] <- unlist(lapply(output_AAFBF_H1, extract_res, 2)) # Bayes factor H0vsH1
        results_H1[, 4] <- unlist(lapply(output_AAFBF_H1, extract_res, 3)) #posterior model probabilities of H0
        
        colnames(results_H1) <- c("BF.10", "PMP.1", "BF.01","PMP.0")
        
        results_H0[, 1] <- unlist(lapply(output_AAFBF_H0, extract_res, 1)) # Bayes factor H1vsH0
        results_H0[, 2] <- unlist(lapply(output_AAFBF_H0, extract_res, 4)) #posterior model probabilities of H1
        results_H0[, 3] <- unlist(lapply(output_AAFBF_H0, extract_res, 2)) # Bayes factor H0vsH1
        results_H0[, 4] <- unlist(lapply(output_AAFBF_H0, extract_res, 3)) #posterior model probabilities of H0
        
        colnames(results_H0) <- c("BF.10", "PMP.1", "BF.01", "PMP.0")
        
        #Evaluation of condition -------------------------------------------
        # Proportions
        prop_BF10 <- length(which(results_H1[, "BF.10"] > BF_thresh)) / ndatasets
        prop_BF01 <- length(which(results_H0[, "BF.01"] > BF_thresh)) / ndatasets
        
        print("Bayes factor check!")
    }
    final_SSD[[b_fract + 1]] <- list(null, hypothesis1)
    final_SSD[[b_fract + 2]] <- BF_thresh
    
    # Final output -----
    print_results(final_SSD)
    if (any(singular_warn > 0)) warning("At least one of the fitted models is singular. For more information about singularity see help('issingular').
                               The number of models that are singular can be found in the output object.")
    invisible(final_SSD)
}
