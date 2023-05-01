######################## DATA GENERATION ##############################

gen_CRT_data <- function(n.datasets = n.datasets, n1 = n1, n2 = n2, var.u0 = var.u0, 
                         var.e = var.e, rho = rho, eff.size = eff.size, hypoth = hypoth, 
                         mean.ctrl = mean.ctrl) {
    ID <- rep(1:n2, each = n1)
    condition <- rep(c(0, 1), each = n1 * n2 / 2)
    # Dummy variables for no intercept model
    Dintervention <- condition
    Dcontrol <- 1 - Dintervention
    mean.interv <- eff.size
    lmer.list <- vector(mode = 'list', length = n.datasets)
    bain.list <- vector(mode = 'list', length = n.datasets)
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
        # bain ------------------------------------------------------------------------
        n.eff <- ((n1 * n2) / (1 + (n1 - 1) * rho.data))/2
        output.bain <- bain(estimates, hypoth, 
                            n = c(n.eff, n.eff), group_parameters = 1, Sigma = cov.list, joint_parameters = 0)
        
        lmer.list[[i]] <- output.lmer
        bain.list[[i]] <- output.bain
        print(i)
    }
    return(output = list(Multilevel = lmer.list,
                         Bayes_factor = bain.list))
}

# Test ---------------------------------------------------------------------------
#a <- gen_CRT_data(10, n1 = 15, n2 = 30, var.u0 = 0.3, var.e = 0.7, rho = 0.3, eff.size = 0.5)

# TODO
# - Change according to the hypothesis.
