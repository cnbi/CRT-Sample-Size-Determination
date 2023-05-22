############################### EVALUATION ####################################


# Evaluation of threshold

eval_thresh <- function(results.H0 = results.H0, results.H1 = results.H1, BF.thresh = BF.thresh, 
                        n.datasets = n.datasets, condition = condition, eta = eta){
    # Proportion
    prop.BF01 <- length(which(results.H0[, 'BF.01'] > BF.thresh)) / n.datasets 
    prop.BF10 <- length(which(results.H1[, 'BF.10'] > BF.thresh)) / n.datasets 
    # Evaluation
    ifelse(prop.BF12 > eta & prop.BF21 > eta, condition <- TRUE, condition <- FALSE)
    #browser()
    #Output
    print(c("prop.BF01: ", prop.BF12, "prop.BF10: ", prop.BF21))
    return(condition)
}

#Test --------------------------------------------------------------------------------