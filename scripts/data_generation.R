######################## DATA GENERATION ##############################

gen_CRT_data <- function(n.datasets = n.datasets, n1 = n1, n2 = n2, var.u0 = var.u0,
                         var.e = var.e, mean.interv, hypoth, b = b) {
    # Create variables ID  of the cluster and condition
    ID <- rep(1:n2, each = n1)
    condition <- rep(c(0, 1), each = n1 * n2 / 2)
    # Dummy variables for no intercept model
    Dintervention <- condition
    Dcontrol <- 1 - Dintervention
    mean.ctrl <- 0
    lmer.list <- vector(mode = 'list', length = n.datasets)
    bain.list <- vector(mode = 'list', length = n.datasets)
    marker <- 0

    # Table for results
    results <- matrix(0, nrow = n.datasets, ncol = 7)

    # simulation of data
    for (i in seq(n.datasets)) {
     # Data generation ----------------------------------------------------------
        set.seed((i + 90) * i)
        u0 <- rnorm(n2, 0, sqrt(var.u0))
        u0 <- rep(u0, each = n1)
        e <- rnorm(n1 * n2, 0, sqrt(var.e))
        resp <- mean.ctrl * Dcontrol + mean.interv * Dintervention + u0 + e

        #Data frame
        data <- cbind(resp, Dintervention, Dcontrol, ID)
        data <- as.data.frame(data)

        # Multilevel analysis ---------------------------------------------------------
        output.lmer <- lmer(resp ~ Dintervention + Dcontrol - 1 + (1 | ID), data = data)
        estimates <- fixef(output.lmer) #means
        cov.intervention <- matrix(vcov(output.lmer)[1], nrow = 1, ncol = 1) #variance-covariance matrix
        cov.control <- matrix(vcov(output.lmer)[4], nrow = 1, ncol = 1)
        cov.list <- list(cov.intervention, cov.control)
        variances <- as.data.frame(VarCorr(output.lmer))
        var.u0.data <- variances[1, 4]
        var.e.data <- variances[2, 4]
        total.var.data <- var.u0.data + var.e.data
        rho.data <- var.u0.data / total.var.data
        ifelse(isSingular(output.lmer), marker <- 1, marker <- 0)

        # bain ------------------------------------------------------------------------
        n.eff <- ((n1 * n2) / (1 + (n1 - 1) * rho.data))/2
        output.bain <- bain(estimates, hypothesis = hypoth,
                            n = c(n.eff, n.eff), group_parameters = 1, Sigma = cov.list,
                            joint_parameters = 0, fraction = b)

        # Results ---------------------------------------------------------------------
        lmer.list[[i]] <- output.lmer
        bain.list[[i]] <- output.bain
        print(c(i, n1, n2)) #change this
        #browser()
        results[i, 1] <- output.bain$estimates[2] # Coefficient control
        results[i, 2] <- output.bain$estimates[1] # Coefficient intervention
        results[i, 3] <- output.bain$BFmatrix[1, 2] # Bayes factor H1vsH2 or H1vsH0
        results[i, 4] <- output.bain$BFmatrix[2, 1] # Bayes factor H2vsH1 or H0vsH1
        results[i, 5] <- output.bain$fit$PMPa[1] #posterior model probabilities of H1.
        results[i, 6] <- output.bain$fit$PMPa[2] #posterior model probabilities of H2 or H0
        results[i, 7] <- marker
    }
    return(output = results)
}




# 
# # Vectorization version ---------------------------------------------------------------------------------
# gen_CRT_data <- function(n.datasets = n.datasets, n1 = n1, n2 = n2, var.u0 = var.u0, 
#                          var.e = var.e, mean.interv, hypoth, b = b) {
#     # Create variables ID  of the cluster and condition
#     ID <- rep(1:n2, each = n1)
#     condition <- rep(c(0, 1), each = n1 * n2 / 2)
#     # Dummy variables for no intercept model
#     Dintervention <- condition
#     Dcontrol <- 1 - Dintervention
#     mean.ctrl <- 0
#     lmer.list <- vector(mode = 'list', length = n.datasets)
#     bain.list <- vector(mode = 'list', length = n.datasets)
#     marker <- 0
#     
#     # Table for results
#     results <- vector(mode = "list", length = n.datasets)
#     data.list <- vector(mode = "list", length = n.datasets)
# 
#     for (i in seq(n.datasets)) {
#         # Data generation ----------------------------------------------------------
#         set.seed((i + 90) * i)
#         u0 <- rnorm(n2, 0, sqrt(var.u0)) 
#         u0 <- rep(u0, each = n1)
#         e <- rnorm(n1 * n2, 0, sqrt(var.e))
#         resp <- mean.ctrl * Dcontrol + mean.interv * Dintervention + u0 + e
#         
#         #Data frame
#         data <- cbind(resp, Dintervention, Dcontrol, ID)
#         data.list[[i]] <- as.data.frame(data)
#     }
#     
#     # Multilevel analysis ---------------------------------------------------------
#     output.lmer <- lapply(data.list, fit_lmer) 
#     estimates <- lapply(output.lmer, fixef)                 # Means
#     cov.intervention <- lapply(output.lmer, varcov, 1)      # Covariance
#     cov.control <- lapply(output.lmer, varcov, 4)
#     cov.list <- Map(list, cov.intervention, cov.control)
#     var.u0.data <- unlist(lapply(output.lmer, get_variances, 1))
#     var.e.data <- unlist(lapply(output.lmer, get_variances, 2))
#     total.var.data <- var.u0.data + var.e.data
#     rho.data <- var.u0.data / total.var.data
#     marker <- lapply(output.lmer, marker_func)
#     
#     browser()
#     #Approximated adjusted fractional Bayes factors 
#     n.eff <- ((n1 * n2) / (1 + (n1 - 1) * rho.data))/2
#     
#     
#     output.bain <- bain(estimates, hypothesis = hypoth, 
#                         n = c(n.eff, n.eff), group_parameters = 1, Sigma = cov.list, 
#                         joint_parameters = 0, fraction = b)
#     
#     # Results ---------------------------------------------------------------------
#     lmer.list[[i]] <- output.lmer
#     bain.list[[i]] <- output.bain
#     print(c(i, n1, n2)) #change this
#     #browser()
#     results[i, 1] <- output.bain$estimates[2] # Coefficient control
#     results[i, 2] <- output.bain$estimates[1] # Coefficient intervention
#     results[i, 3] <- output.bain$BFmatrix[1, 2] # Bayes factor H1vsH2 or H1vsH0
#     results[i, 4] <- output.bain$BFmatrix[2, 1] # Bayes factor H2vsH1 or H0vsH1
#     results[i, 5] <- output.bain$fit$PMPa[1] #posterior model probabilities of H1.
#     results[i, 6] <- output.bain$fit$PMPa[2] #posterior model probabilities of H2 or H0
#     results[i, 7] <- marker
#     
#     return(output = results)
# }


# Test ---------------------------------------------------------------------------
# a <- gen_CRT_data(10, n1 = 15, n2 = 30, var.u0 = 0.3, var.e = 0.7, mean.interv = 0.4, hypoth = hypoth, b = 1)
