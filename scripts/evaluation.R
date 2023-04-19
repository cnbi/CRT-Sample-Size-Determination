############################### EVALUATION ####################################


# Evaluation of threshold

eval_thresh <- function(results.H0 = results.H0, results.H1 = results.H1, BF.thresh = BF.thresh, n.datasets = n.datasets, condition = condition, eta = eta){
    # Proportion
    prop.BF12 <- length(which(results.H0[, 'BF.12'] > BF.thresh)) / n.datasets #Or 01 instead? Proportion of BF12 when H0 = TRUE
    prop.BF21 <- length(which(results.H1[, 'BF.21'] > BF.thresh)) / n.datasets #Or 10? Proportion of BF21 when H1 = TRUE
    # Evaluation
    ifelse(prop.BF12 > eta & prop.BF21 > eta, condition <- TRUE, condition <- FALSE)
    #browser()
    #Output
    print(c("prop.BF12: ", prop.BF12, "prop.BF21: ", prop.BF21))
    return(condition)
}

#Test --------------------------------------------------------------------------------