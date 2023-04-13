############################### EVALUATION ####################################


# Evaluation of threshold

eval.thresh <- function(results, BF.thresh = BF.thresh, n.datasets = n.datasets){
    # Proportion
    prop.BF12 <- length(which(results[, 'BF.12'] > BF.thresh)) / n.datasets #Or 01 instead?
    prop.BF21 <- length(which(results[, 'BF.21'] > BF.thresh)) / n.datasets #Or 10?
    
    #Output
    return(list(prop.BF12 = prop.BF12,
                prop.BF21 = prop.BF21))
}
