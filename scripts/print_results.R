################# PRINT RESULTS #############################

print_results <- function(object_result) {
    title <- "Final sample size"
    cat(paste("\n", title, "\n", sep = ""))
    row <- paste(rep("=", nchar(title)), collapse = "")
    cat(row, "\n")
    if (object_result[[length(object_result)]] == "Informative") {   # Print for informative hypotheses
        cat("Hypotheses:", "\n")
        cat("    H1:", object_result[[5]][[1]], "\n")
        cat("    H2:", object_result[[5]][[2]], "\n")
        cat("Using cluster size = ", object_result$n1, " and number of clusters = ", object_result$n2, "\n")
        cat("P (BF.12 > ", object_result[[6]], " | H1) = ", object_result$Eta, "\n")
    } else { # Print for null vs informative
        n_object <- length(object_result)
        results_matrix <- matrix(NA, nrow = length(object_result) - 2, ncol = 5)
        results_matrix[, 1] <- seq(length(object_result) - 2)
        results_matrix[, 2] <- unlist(lapply(object_result, `[[`, 1))
        results_matrix[, 3] <- unlist(lapply(object_result, `[[`, 2))
        results_matrix[, 4] <- unlist(lapply(object_result, `[[`, 3))
        results_matrix[, 5] <- unlist(lapply(object_result, `[[`, 4))
        rownames(results_matrix) <- c("b", "n1", "n2", paste("P(BF.01 >", object_result[[n_object]], "| H0)", sep = " "), 
                                      paste("P(BF.10 >", object_result[[n_object]], "| H1)", sep = " "))
        
        cat("Hypotheses:", "\n")
        cat("    H0:", object_result[[4]][[1]], "\n")
        cat("    H1:", object_result[[4]][[2]], "\n")
        
        cat("******************************", "\n")
        print(format(results_matrix, justify = "centre"))
        cat("\n", "******************************", "\n")
        cat("n1: Cluster size", "\n")
        cat("n2: Number of clusters")
    }
}