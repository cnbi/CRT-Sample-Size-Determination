######################## DATA GENERATION ##############################

# Vectorization version ---------------------------------------------------------------------------------
# gen_CRT_data <- function(ndatasets = ndatasets, n1 = n1, n2 = n2, var_u0 = var_u0,
#                          var_e = var_e, mean_interv, b = b, type, batch_size) {
#     # Create variables id  of the cluster and condition
#     id <- rep(1:n2, each = n1)
#     condition <- rep(c(0, 1), each = n1 * n2 / 2)
#     # Dummy variables for no intercept model
#     intervention <- condition
#     control <- 1 - intervention
#     mean_control <- 0
#     marker <- 0
#     
#     # Tables for results
#     results <- matrix(NA, ndatasets, 4)
#     output_lmer <- vector(mode = "list", length = ndatasets)
#     data.list <- vector(mode = "list", length = ndatasets)
#     
#     # Data generation ----------------------------------------------------------
#     for (iter in seq(ndatasets)) {
#         set.seed((iter + 90) * iter)
#         u0 <- rnorm(n2, 0, sqrt(var_u0))
#         u0 <- rep(u0, each = n1)
#         e <- rnorm(n1 * n2, 0, sqrt(var_e))
#         resp <- mean_control * control + mean_interv * intervention + u0 + e
#         #Data frame
#         data.list[[iter]] <- cbind(resp, intervention, control, id)
#     }
#     data.list <- lapply(data.list, as.data.frame) # Necessary to use lmer
# 
#     # Multilevel analysis ---------------------------------------------------------
#     # Batches
#     batch_size <- batch_size
#     ifelse((ndatasets/batch_size) %% 1 == 0, num_batches <- ndatasets/batch_size, 
#            num_batches <- (ndatasets/batch_size) + 1)
#     for (batch in seq(num_batches)) {
#         #Indexes
#         start_index <- (batch_size * (batch - 1)) + 1
#         end_index <- min(batch * batch_size, ndatasets)
#         #Multilevel fitting
#         output_lmer[start_index:end_index] <- lapply(data.list[start_index:end_index], fit_lmer)
#     }
#     marker <- lapply(output_lmer, marker_func) # Mark singularity
#     
#     estimates <- lapply(output_lmer, fixef)                 # Means
#     cov_intervention <- lapply(output_lmer, varcov, 1)      # Covariance
#     cov_control <- lapply(output_lmer, varcov, 4)
#     cov_list <- Map(list, cov_intervention, cov_control)
#     var_u0_data <- unlist(lapply(output_lmer, get_variances, 1))
#     var_e_data <- unlist(lapply(output_lmer, get_variances, 2))
#     total_var_data <- var_u0_data + var_e_data
#     rho_data <- var_u0_data / total_var_data
#     print("Multilevel check")
#   
#     #Approximated adjusted fractional Bayes factors------------------------------
#     n_eff <- ((n1 * n2) / (1 + (n1 - 1) * rho_data)) / 2
#     output_AAFBF <- Map(calc_aafbf, type, estimates, cov_list, list(b), n_eff)
# 
#     # Results ---------------------------------------------------------------------
#     results[, 1] <- unlist(lapply(output_AAFBF, extract_res, 1)) # Bayes factor H1vsH2 or H1vsH0
#     results[, 2] <- unlist(lapply(output_AAFBF, extract_res, 2)) # Bayes factor H2vsH1 or H0vsH1
#     results[, 3] <- unlist(lapply(output_AAFBF, extract_res, 3)) #posterior model probabilities of H2 or H0
#     results[, 4] <- unlist(lapply(output_AAFBF, extract_res, 4)) #posterior model probabilities of H1
#     rm(id, condition, intervention, control, mean_control, data.list, output_lmer, u0,
#        e, resp, estimates, batch_size, cov_intervention, cov_control, cov_list,  var_u0_data, 
#        var_e_data, total_var_data, rho_data, n_eff, output_AAFBF)
#     
#     return(output = results)
# }


# Test ---------------------------------------------------------------------------
#a <- gen_CRT_data(10, n1 = 15, n2 = 30, var_u0 = 0.3, var_e = 0.7,
#                  mean_interv = 0.4, b = 2, type = "equality")


# New version -------------------------------------------------------------------
# Vectorization version ---------------------------------------------------------------------------------
gen_CRT_data <- function(ndatasets = ndatasets, n1 = n1, n2 = n2, var_u0 = var_u0,
                         var_e = var_e, mean_interv, batch_size) {
  # Create variables id  of the cluster and condition
  id <- rep(1:n2, each = n1)
  condition <- rep(c(0, 1), each = n1 * n2 / 2)
  # Dummy variables for no intercept model
  intervention <- condition
  control <- 1 - intervention
  mean_control <- 0
  marker <- 0
  
  # Tables for results
  results <- matrix(NA, ndatasets, 4)
  output_lmer <- vector(mode = "list", length = ndatasets)
  data.list <- vector(mode = "list", length = ndatasets)
  
  # Data generation ----------------------------------------------------------
  for (iter in seq(ndatasets)) {
    set.seed((iter + 90) * iter)
    u0 <- rnorm(n2, 0, sqrt(var_u0))
    u0 <- rep(u0, each = n1)
    e <- rnorm(n1 * n2, 0, sqrt(var_e))
    resp <- mean_control * control + mean_interv * intervention + u0 + e
    #Data frame
    data.list[[iter]] <- cbind(resp, intervention, control, id)
  }
  data.list <- lapply(data.list, as.data.frame) # Necessary to use lmer
  
  # Multilevel analysis ---------------------------------------------------------
  # Batches
  batch_size <- batch_size
  ifelse((ndatasets/batch_size) %% 1 == 0, num_batches <- ndatasets/batch_size, 
         num_batches <- (ndatasets/batch_size) + 1)
  for (batch in seq(num_batches)) {
    #Indexes
    start_index <- (batch_size * (batch - 1)) + 1
    end_index <- min(batch * batch_size, ndatasets)
    #Multilevel fitting
    output_lmer[start_index:end_index] <- lapply(data.list[start_index:end_index], fit_lmer)
  }
  marker <- lapply(output_lmer, marker_func) # Mark singularity
  singular_datasets <- Reduce("+", marker)
  
  estimates <- lapply(output_lmer, fixef)                 # Means
  cov_intervention <- lapply(output_lmer, varcov, 1)      # Covariance
  cov_control <- lapply(output_lmer, varcov, 4)
  cov_list <- Map(list, cov_intervention, cov_control)
  var_u0_data <- unlist(lapply(output_lmer, get_variances, 1))
  var_e_data <- unlist(lapply(output_lmer, get_variances, 2))
  total_var_data <- var_u0_data + var_e_data
  rho_data <- var_u0_data / total_var_data
  print("Multilevel check")

  rm(id, condition, intervention, control, mean_control, data.list, output_lmer, u0,
     e, resp, batch_size, cov_intervention, cov_control, var_u0_data, 
     var_e_data, total_var_data)
  
  return(output <- list("rho_data" = rho_data,
                        "estimates" = estimates,
                        "cov_list" = cov_list,
                        "singularity" = singular_datasets))
}
