############################ EXTRACT RESULTS #############################

extract_results <- function(n.datasets = n.datasets, data, names.table = names.table){
    results <- matrix(nrow = n.datasets, ncol = 11)
    # Extract Bayes Factors
    for (i in seq(n.datasets)) {
        results[i, 1] <- data[[2]][[i]]$estimates[2] # Coefficient control
        results[i, 2] <- data[[2]][[i]]$estimates[1] # Coefficient intervention
        results[i, 3] <- data[[2]][[i]]$BFmatrix[1, 2] # Bayes factor H1vsH2
        results[i, 4] <- data[[2]][[i]]$BFmatrix[2, 1] # Bayes factor H2vsH1
        results[i, 5] <- data[[2]][[i]]$fit[1, 7] # Bayes factor of H1 versus its complement.
        results[i, 6] <- data[[2]][[i]]$fit[2, 7] # Bayes factor of H2 versus its complement.
        results[i, 7] <- data[[2]][[i]]$fit$PMPa[1] #posterior model probabilities of H1.
        results[i, 8] <- data[[2]][[i]]$fit$PMPa[2] #posterior model probabilities of H2.
        results[i, 9] <- data[[2]][[i]]$fit$PMPc[1] #posterior model probabilities of H1 + its complement.
        results[i, 10] <- data[[2]][[i]]$fit$PMPc[2] #posterior model probabilities of H2 + its complement.
        results[i, 11] <- data[[2]][[i]]$fit$PMPc[4] #posterior model probabilities of Hc
    }
    colnames(results) <- names.table
    # Output
    return(results)
}

# Test -----------------------------------------------------------------------
#a.results <- extract_results(n.datasets = 10, data = a)
