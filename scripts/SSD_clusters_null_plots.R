############# Sample Size Determination for Cluster Randomized Trials #############



# FOR PLOTS
# With this function is possibl to generate the data necessary to make plots.


SSD_crt_null_plots <- function(eff.size, n1 = 15, n2 = 30, n.datasets = 1000, rho, BF.thresh, 
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
    best_result <- FALSE
    
    maxim <- 200
    
    #Hypotheses -----
    hypothesis1 <- "Dintervention>Dcontrol"
    null <- "Dintervention=Dcontrol"
    hypoth <- paste(null, ";", hypothesis1)
    final_SSD <- vector(mode = "list", length = 3)
    
    SSD_object <- vector(mode = "list", maxim - n2)
    
    # Simulation and evaluation of condition  -----
    for (b in seq(b.fract)) {
        while (n2 < maxim + 1) {
            # If H1 is true
            data.H1 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e, 
                                                  mean.interv = eff.size, hypoth, b))
            colnames(data.H1) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
            
            # If H0 is true
            data.H0 <- do.call(gen_CRT_data, list(n.datasets, n1, n2, var.u0, var.e,
                                                  mean.interv = eff.size0, hypoth, b))
            colnames(data.H0) <- c('Dcontrol', 'Dintervention', 'BF.01', 'BF.10', 'PMP.0', 'PMP.1')
            
            #Evaluation of condition ----
            # Proportion
            prop.BF01 <- length(which(data.H0[, 'BF.01'] > BF.thresh)) / n.datasets 
            prop.BF10 <- length(which(data.H1[, 'BF.10'] > BF.thresh)) / n.datasets 
            # Evaluation
            ifelse(prop.BF01 > eta & prop.BF10 > eta, condition <- TRUE, condition <- FALSE)
            
            # Output ----
            
            print(c("Using cluster size: ", n1, "and number of clusters: ", n2,
                    "prop.BF01: ", prop.BF01, "prop.BF10: ", prop.BF10))
            if (fixed == 'n1') {
                n2 = n2 + 2
            } else if (fixed == 'n2') {
                n1 = n1 + 1
            }
            
            iterations <- iterations + 1
            SSD_object[[iterations]] <- list("n1" = n1,
                                             "n2" = n2,
                                             "Proportion.BF01" = prop.BF01,
                                             "Proportion.BF10" = prop.BF10,
                                             "b.frac" = b,
                                             "data.H0" = data.H0,
                                             "data.H1" = data.H1)
            
        }
        print("Finished!")
        
    }
    
    # Final output -----
    
    return(SSD_object)
}


