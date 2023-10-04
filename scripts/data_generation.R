######################## DATA GENERATION ##############################

# gen_CRT_data <- function(ndatasets = ndatasets, n1 = n1, n2 = n2, var_u0 = var_u0,
#                          var_e = var_e, mean_interv, hypoth, b = b) {
#     # Create variables id  of the cluster and condition
#     id <- rep(1:n2, each = n1)
#     condition <- rep(c(0, 1), each = n1 * n2 / 2)
#     # Dummy variables for no intercept model
#     intervention <- condition
#     control <- 1 - intervention
#     mean_control <- 0
#     lmer_list <- vector(mode = "list", length = ndatasets)
#     bain_list <- vector(mode = "list", length = ndatasets)
#     marker <- 0
# 
#     # Table for results
#     results <- matrix(0, nrow = ndatasets, ncol = 7)
# 
#     # simulation of data
#     for (i in seq(ndatasets)) {
#      # Data generation --------------------------------------------------------
#         set.seed((i + 90) * i)
#         u0 <- rnorm(n2, 0, sqrt(var_u0))
#         u0 <- rep(u0, each = n1)
#         e <- rnorm(n1 * n2, 0, sqrt(var_e))
#         resp <- mean_control * control + mean_interv * intervention + u0 + e
# 
#         #Data frame
#         data <- cbind(resp, intervention, control, id)
#         data <- as.data.frame(data)
# 
#         # Multilevel analysis --------------------------------------------------
#         output_lmer <- lmer(resp ~ intervention + control - 1 + (1|id), data = data)
#         estimates <- fixef(output_lmer)
#         cov_intervention <- matrix(vcov(output_lmer)[1], nrow = 1, ncol = 1) #variance-covariance matrix
#         cov_control <- matrix(vcov(output_lmer)[4], nrow = 1, ncol = 1)
#         cov_list <- list(cov_intervention, cov_control)
#         variances <- as.data.frame(VarCorr(output_lmer))
#         var_u0_data <- variances[1, 4]
#         var_e_data <- variances[2, 4]
#         total_var_data <- var_u0_data + var_e_data
#         rho_data <- var_u0_data / total_var_data
#         ifelse(isSingular(output_lmer), marker <- 1, marker <- 0)
# 
#         # bain -----------------------------------------------------------------
#         n_eff <- ((n1 * n2) / (1 + (n1 - 1) * rho_data)) / 2
#         output_bain <- bain(estimates, hypothesis = hypoth,
#                             n = c(n_eff, n_eff), group_parameters = 1, Sigma = cov_list,
#                             joint_parameters = 0, fraction = b)
# 
#         # Results --------------------------------------------------------------
#         lmer_list[[i]] <- output_lmer
#         bain_list[[i]] <- output_bain
#         print(c(i, n1, n2)) #change this
#         #browser()
#         results[i, 1] <- output_bain$estimates[2] # Coefficient control
#         results[i, 2] <- output_bain$estimates[1] # Coefficient intervention
#         results[i, 3] <- output_bain$BFmatrix[1, 2] # Bayes factor H1vsH2 or H1vsH0
#         results[i, 4] <- output_bain$BFmatrix[2, 1] # Bayes factor H2vsH1 or H0vsH1
#         results[i, 5] <- output_bain$fit$PMPa[1] #posterior model probabilities of H1.
#         results[i, 6] <- output_bain$fit$PMPa[2] #posterior model probabilities of H2 or H0
#         results[i, 7] <- marker
#     }
#     return(output = results)
# }



# Vectorization version ---------------------------------------------------------------------------------
gen_CRT_data <- function(ndatasets = ndatasets, n1 = n1, n2 = n2, var_u0 = var_u0,
                         var_e = var_e, mean_interv, b = b, type) {
    # Create variables id  of the cluster and condition
    id <- rep(1:n2, each = n1)
    condition <- rep(c(0, 1), each = n1 * n2 / 2)
    # Dummy variables for no intercept model
    intervention <- condition
    control <- 1 - intervention
    mean_control <- 0
    lmer_list <- vector(mode = "list", length = ndatasets)
    bayes_list <- vector(mode = "list", length = ndatasets)
    marker <- 0
    # without_sing <- FALSE
    
    # Table for results
    results <- matrix(NA, ndatasets, 4)
    data.list <- vector(mode = "list", length = ndatasets)
    output_lmer <- vector(mode = "list", length = ndatasets)
    
    for (iter in seq(ndatasets)) {
        # Data generation ----------------------------------------------------------
        set.seed((iter + 90) * iter)
        u0 <- rnorm(n2, 0, sqrt(var_u0))
        u0 <- rep(u0, each = n1)
        e <- rnorm(n1 * n2, 0, sqrt(var_e))
        resp <- mean_control * control + mean_interv * intervention + u0 + e
        
        #Data frame
        data <- cbind(resp, intervention, control, id)
        data.list[[iter]] <- as.data.frame(data)
    }
    
    # Multilevel analysis ---------------------------------------------------------
    output_lmer <- lapply(data.list, fit_lmer)
    marker <- lapply(output_lmer, marker_func) # Mark singularity
    # ifelse(sum(unlist(marker)) > 0, without_sing <- FALSE, without_sing <- TRUE)
    
    estimates <- lapply(output_lmer, fixef)                 # Means
    cov_intervention <- lapply(output_lmer, varcov, 1)      # Covariance
    cov_control <- lapply(output_lmer, varcov, 4)
    cov_list <- Map(list, cov_intervention, cov_control)
    var_u0_data <- unlist(lapply(output_lmer, get_variances, 1))
    var_e_data <- unlist(lapply(output_lmer, get_variances, 2))
    total_var_data <- var_u0_data + var_e_data
    rho_data <- var_u0_data / total_var_data
    #Approximated adjusted fractional Bayes factors
    #browser()
    n_eff <- ((n1 * n2) / (1 + (n1 - 1) * rho_data)) / 2
    output_AAFBF <- Map(calc_aafbf, type, estimates, cov_list, list(b), n_eff)
    # output_bain <- bain(estimates, hypothesis = hypoth,
    #                      n = c(n_eff, n_eff), group_parameters = 1, Sigma = cov_list,
    #                      joint_parameters = 0, fraction = b)
    
    # Results ---------------------------------------------------------------------
    results[, 1] <- unlist(lapply(output_AAFBF, extract_res, 1)) # Bayes factor H1vsH2 or H1vsH0
    results[, 2] <- unlist(lapply(output_AAFBF, extract_res, 2)) # Bayes factor H2vsH1 or H0vsH1
    results[, 3] <- unlist(lapply(output_AAFBF, extract_res, 3)) #posterior model probabilities of H2 or H0
    results[, 4] <- unlist(lapply(output_AAFBF, extract_res, 4)) #posterior model probabilities of H1.
    
    
    return(output = results)
}


# Test ---------------------------------------------------------------------------
#a <- gen_CRT_data(10, n1 = 15, n2 = 30, var_u0 = 0.3, var_e = 0.7, mean_interv = 0.4, hypoth = hypoth, b = 1, type = "equality")
