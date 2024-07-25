# Function to compare b

# Check whether the sample size increases with the increase of fraction b
# Condition: Increases with fraction b.
# 

check_b <- function(nrow_design, dataframe, design_matrix){
    # Filter data set with final results for only one condition
    filtered_df <- dataframe[which(dataframe$rho == design_matrix[nrow_design, 1] &
                                       dataframe$eff_size == design_matrix[nrow_design, 2] &
                                       dataframe$BF_threshold == design_matrix[nrow_design, 3]), ]
    if (dataframe$fixed[1] == "n2") {
        filtered_df <- filtered_df[which(filtered_df$n2 == design_matrix[nrow_design, 4]), ]
    } else {
        filtered_df <- filtered_df[which(filtered_df$n1 == design_matrix[nrow_design, 4]), ]
    }
    
    # Check that sample size for b1 < b2 < b3
    if (dataframe$fixed[1] == "n2") {
        ifelse(filtered_df[which(filtered_df$b == 1), 16] < filtered_df[which(filtered_df$b == 2), 16] && filtered_df[which(filtered_df$b == 2), 16] < filtered_df[which(filtered_df$b == 3), 16], increase <- TRUE, increase <- FALSE) 
    } else {
        ifelse(filtered_df[which(filtered_df$b == 1), 15] < filtered_df[which(filtered_df$b == 2), 15] && filtered_df[which(filtered_df$b == 2), 15] < filtered_df[which(filtered_df$b == 3), 15], increase <- TRUE, increase <- FALSE)
    }
    
    # Check that sample size for b1 = b2 = b3
    if (dataframe$fixed[1] == "n2") {
        ifelse(filtered_df[which(filtered_df$b == 1), 16] == filtered_df[which(filtered_df$b == 2), 16] && filtered_df[which(filtered_df$b == 2), 16] == filtered_df[which(filtered_df$b == 3), 16], maintain <- TRUE, maintain <- FALSE) 
    } else {
        ifelse(filtered_df[which(filtered_df$b == 1), 15] == filtered_df[which(filtered_df$b == 2), 15] && filtered_df[which(filtered_df$b == 2), 15] == filtered_df[which(filtered_df$b == 3), 15], maintain <- TRUE, maintain <- FALSE)
    }
    
    # Check that sample size for b1 > b2 > b3
    if (dataframe$fixed[1] == "n2") {
        ifelse(filtered_df[which(filtered_df$b == 1), 16] > filtered_df[which(filtered_df$b == 2), 16] && filtered_df[which(filtered_df$b == 2), 16] > filtered_df[which(filtered_df$b == 3), 16], decrease <- TRUE, decrease <- FALSE) 
    } else {
        ifelse(filtered_df[which(filtered_df$b == 1), 15] > filtered_df[which(filtered_df$b == 2), 15] && filtered_df[which(filtered_df$b == 2), 15] > filtered_df[which(filtered_df$b == 3), 15], decrease <- TRUE, decrease <- FALSE)
    }
    
    # Return result
    return(c(increase, maintain, decrease))
}

# Finding n1
result_b_N1 <- matrix(NA, nrow = nrow_designN1, ncol = 3)
colnames(result_b_N1) <- c("increased", "maintained", "decreased")
for (ROW in 1:nrow_designN1) {
    # Evaluate
    result_b_N1[ROW, ] <- check_b(ROW, dataframe = final_results_FindN1, design_matrix = design_matrixN1)
}
result_b_N1 <- cbind(design_matrixN1, result_b_N1)

# Finding n2
result_b_N2 <- matrix(NA, nrow = nrow_designN1, ncol = 3)
colnames(result_b_N2) <- c("increased", "maintained", "decreased")
for (ROW in 1:nrow_designN2) {
    # Evaluate
    result_b_N2[ROW, ] <- check_b(ROW, dataframe = final_results_FindN2, design_matrix = design_matrixN2)
}
result_b_N2 <- cbind(design_matrixN2, result_b_N2)

decreases_N2 <- result_b_N2[which(result_b_N2$decreased == TRUE), ]
others_N2 <- result_b_N2[which((result_b_N2$increased == FALSE) & (result_b_N2$maintained == FALSE) & (result_b_N2$decreased == FALSE)), ]
frequency_tables <- lapply(decreases_N2[names(result_b_N2)], table)

# Check frequency tables
print(frequency_tables)

decreases_N1 <- result_b_N1[which(result_b_N1$decreased == TRUE), ]
others_N1 <- result_b_N1[which((result_b_N1$increased == FALSE) & (result_b_N1$maintained == FALSE) & (result_b_N1$decreased == FALSE)), ]
frequency_tables <- lapply(others_N1[names(result_b_N1)], table)
