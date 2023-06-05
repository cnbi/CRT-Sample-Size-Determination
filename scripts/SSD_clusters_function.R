############# Sample Size Determination for Cluster Randomized Trials #############

#'@title Sample Size Determination for Cluster Randomized Trials
#'@description 
#'@arguments
## eff.size: Effect size
## n.datasets: Number of data sets that the user wants to generate to determine the sample size.
## var.u0: Between cluster variance.
## BF.thresh: Threshold to choose.
## eta: Posterior probability that is expected. 1-eta = Bayesian error probability.
## plots: Logical argument. If TRUE the function print the plots.
## hypothesis: The hypothesis that is going to be tested.


SSD_crt_null <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh, 
                    eta = 0.8, fixed = 'n2', b.fract = 3) {
    # Libraries ----
    library(lme4)
    library(bain)
    library(dplyr)
    
    #Functions ----
    source('data_generation.R')
    # source('print.null')
    
    # Starting values -----
    total.var <- 1
    var.u0 <- rho * total.var
    var.e <- total.var - var.u0
    iterations <- 1
    eff.size0 <- 0
    condition <- FALSE #condition fulfillment indicator
    
    # Binary search start ----
    # if (fixed == 'n1') {
    #     l <- n2
    # } else if (fixed == 'n2') {
    #     l <- n1
    # }
    # h <- 1000
    
    #Hypotheses -----
    hypothesis1 <- "Dintervention>Dcontrol"
    null <- "Dintervention=Dcontrol"
    hypoth <- paste(null, ";", hypothesis1)
    final_SSD <- vector(mode = "list", length = 3)
    
    # Simulation and evaluation of condition  -----
    for (b in seq(b.fract)) {
        while (condition == FALSE) {
            # If H1 is true
            data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, 
                                                  eff.size, hypoth))
            colnames(data.H1) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
            

            # If H0 is true
            data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e,
                                                  mean.interv = eff.size0, hypoth))
            colnames(data.H0) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
            
            #Evaluation of condition ----
            # Proportion
            prop.BF01 <- length(which(data.H0[, 'BF.01'] > BF.thresh)) / n.datasets 
            prop.BF10 <- length(which(data.H1[, 'BF.10'] > BF.thresh)) / n.datasets 
            # Evaluation
            ifelse(prop.BF01 > eta & prop.BF10 > eta, condition <- TRUE, condition <- FALSE)
            #browser()
            # Binary search algorithm -----
            # if (condition == FALSE) {
            #     print("Using cluster size:", n1, "and number of clusters:", n2, 
            #           "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10, sep = " ")
            #     if (fixed == 'n1') {
            #         l <- n2     #lower bound
            #         h <- h   #higher bound
            #         n2 <- round((l + h) / 2)    #point in the middle
            #     } else if (fixed == 'n2') {
            #         l <- n1     #lower bound
            #         h <- h    #higher bound
            #         n1 <- round((l + h) / 2)    #point in the middle
            #     }
            # } else if (condition == TRUE) {
            #     print("Using cluster size:", n1, "and number of clusters:", n2, 
            #           "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10, sep = " ")
            #     if (fixed == 'n1') {
            #         l <- l    #lower bound
            #         h <- n2    #higher bound
            #         n2 <- round((l + h) / 2)    #point in the middle
            #     } else if (fixed == 'n2') {
            #         l <- l     #lower bound
            #         h <- n1    #higher bound
            #         n1 <- round((l + h) / 2)    #point in the middle
            #     }
            # }
            # Output ----
            if (condition == FALSE) {
                print(c("Using cluster size: ", n1, "and number of clusters: ", n2, 
                      "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10))
                if (fixed == 'n1') {
                    n2 = n2 + 2
                } else if (fixed == 'n2') {
                    n1 = n1 + 1
                }
            }
            
            iterations <- iterations + 1
            if (n2 == 1000) {
                break
            }
            
        }
        
        SSD_object <- list("n1" = n1,
                           "n2" = n2,
                           "Proportion.BF01" = prop.BF01,
                           "Proportion.BF10" = prop.BF10,
                           "b.frac" = b,
                           "data.H0" = data.H0,
                           "data.H1" = data.H1)
        final_SSD[[b]] <- SSD_object
    }
    
    # Final output -----
    for (b in seq(b.fract)) {
        title <- "Final sample size"
        cat(paste("\n", title, "\n", sep = ""))
        row <- paste(rep("=", nchar(title)), collapse = "")
        cat(row, "\n")
        cat("Hypotheses:", "\n")
        cat("H0:", null, "\n")
        cat("H1:", hypothesis1, "\n")
        cat("Using b fraction = ", final_SSD[[b]]$b.frac, 
            ", cluster size = ", final_SSD[[b]]$n1, 
            ", and number of clusters = ", final_SSD[[b]]$n2, "\n")
        cat("P (BF.01 >", BF.thresh, " | H0) = ", final_SSD[[b]]$Proportion.BF01, "\n")
        cat("P (BF.10 >", BF.thresh, " | H1) = ", final_SSD[[b]]$Proportion.BF10, "\n")
    }
    
    return(final_SSD)
}



#TODO:
# - Check style with lintR
# - Think a better name for the condition object.
# - Run a simulation to know see the performance under various conditions.


# Warnings -------------------------------------------------------------------
# Check that n2 is even.
# Check that n1.fixed and n2.fixed are not TRUE both.

# Test -------------------------------------------------------------------------
start.time <- Sys.time()
nulla <- SSD_crt_null(eff.size = 0.5, n.datasets = 10, rho = 0.1, BF.thresh = 3, fixed = "n1", 
             b.fract = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE) #singularity
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.01, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
# This is weird, it seems like this needs less number of clusters.

start.time <- Sys.time()
SSD_crt(eff.size = 0.4, n.datasets = 15, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
try <- SSD_crt_null(eff.size = 0.2, n.datasets = 20, rho = 0.1, BF.thresh = 3,fixed = "n1", 
        b.fract = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.2, n.datasets = 20, rho = 0.05, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
SSD_crt(eff.size = 0.8, n.datasets = 20, rho = 0.1, BF.thresh = 3, hypothesis = "interv.bigger",
        n1.fixed = TRUE, n2.fixed = FALSE)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

