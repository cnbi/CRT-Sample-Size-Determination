###########################TEST INFORMATIVE HYPOTHESIS ############################

# Design matrix
n1 <- c(5, 10, 20, 40)
n2 <- c(30, 60 ,90)
rho <- c(0.25, 0.05, 0.1) #Intraclass correlation
eff.size <- c(0.2, 0.5, 0.8)
bf.thresh <- c(1, 3, 5)
fix <- c("n1", "n2")
design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh)

nrow.design <- nrow(design.matrix)


for (Row in seq(nrow.design)) {
    # function
    ssd_results <- SSD_crt_inform(eff.size = design.matrix[Row, 4], 
                                  n1 = design.matrix[Row, 1],
                                  n2 = design.matrix[Row, 2],
                                  n.datasets = 1000,
                                  rho = design.matrix[Row, 3],
                                  BF.thresh = design.matrix[Row, 5],
                                  eta = 0.8, 
                                  fixed = design.matrix[Row, 6])
    # Save results
    save(ssd_results, file = paste("ResultRow", Row, ".Rdata", sep = ""))

}

# Plots


# - Add plots.This could be a function.