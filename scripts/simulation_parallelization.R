############### FINAL SIMULATION FIRST PROJECT SAVING EACH INCREASE############

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots
library(parallel)
library(foreach)
library(doParallel)
library(pbapply)

# Functions --------------------------------------------------------------------
source("SSD_clusters_function_pl.R")

# Creation of folder for results
path <- "~/"
results_folder <- "results"
dir.create(results_folder)

# Design matrix-----------------------------------------------------------------

rho <- c(0.025, 0.05, 0.1)
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)
b_fract <- 3
# fixed <- c("n1", "n2")
# Fixed n1
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(rho, eff_size, BF_threshold, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
timesN2 <- matrix(NA, nrow = nrow_designN2, ncol = 1)
colnames(design_matrixN2) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

#fixed n2
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
nrow_designN1 <- nrow(design_matrixN1)
timesN1 <- matrix(NA, nrow = nrow_designN1, ncol = 1)
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

# ~N1: To determine n1
# ~N2: To determine n2

# Parallelization -----------------------------------------------------------
## Fixed n1
sim_find_n2 <- function(Row, design_matrixN2, results_folder, b_fract){
    # Start time
    start <- Sys.time()
    # Actual simulation
    ssd_results <- SSD_crt_null_plots(eff_size = design_matrixN2[Row, 2],
                                n1 = design_matrixN2[Row, 4],
                                ndatasets = 5, rho = design_matrixN2[Row, 1],
                                BF_thresh = design_matrixN2[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN2[Row, 5]), b_fract = b_fract)
    # Save results
    timesN2_row <- Sys.time() - start
    
    file_name <- file.path(paste0(results_folder, "/ResultsN2Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    return(timesN2_row)
}
registerDoParallel(2) #Change this
timesN2 <- foreach(Row = seq_len(nrow_designN2), .combine = "c") %dopar% {
    sim_find_n2(Row, design_matrixN2, results_folder, 3)
}
stopImplicitCluster()
    
    
# for (Row in seq(nrow_designN2)) {
#     # Start time
#     start <- Sys.time()
#     # Actual simulation
#     ssd_results <- SSD_crt_null(eff_size = design_matrixN2[Row, 2],
#                                 n1 = design_matrixN2[Row, 4],
#                                 ndatasets = 5000, rho = design_matrixN2[Row, 1],
#                                 BF_thresh = design_matrixN2[Row, 3], eta = 0.8,
#                                 fixed = design_matrixN2[Row, 5], b_fract = b_fract)
#     # Save results
#     timesN2[Row] <- Sys.time() - start
#     file_name <- paste0(results_folder, "/ResultsN2Row", Row, ".RDS")
#     saveRDS(ssd_results, file = file_name)
# }
# timesN2 <- as.data.frame(cbind(design_matrixN2, timesN2))
# saveRDS(timesN2, file = "results/timesN2.RDS")


## Fixed n2
sim_find_n1 <- function(Row, design_matrixN1, results_folder, b_fract){
    # Start time
    start <- Sys.time()
    # Actual simulation
    ssd_results <- SSD_crt_null_plots(eff_size = design_matrixN1[Row, 2],
                                n2 = design_matrixN1[Row, 4],
                                ndatasets = 5000, rho = design_matrixN1[Row, 1],
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = design_matrixN1[Row, 5], b_fract = b_fract)
    # Save results
    timesN1_row <- Sys.time() - start
    file_name <- file.path(paste0(results_folder, "/ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    return(timesN1_row)
}
registerDoParallel(10)
timesN2 <- foreach(Row = seq_len(nrow_designN1), .combine = "c") %dopar% {
    sim_find_n1(Row, design_matrixN1, results_folder, 3)
}
stopImplicitCluster()
timesN1 <- as.data.frame(cbind(design_matrixN1, timesN1))
saveRDS(timesN1, file = "results/timesN1.RDS")
# for (Row in seq(nrow_designN1)) {
#     # Start time
#     start <- Sys.time()
#     # Actual simulation
#     ssd_results <- SSD_crt_null(eff_size = design_matrixN1[Row, 2],
#                                 n2 = design_matrixN1[Row, 4],
#                                 ndatasets = 5000, rho = design_matrixN1[Row, 1],
#                                 BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
#                                 fixed = design_matrixN1[Row, 5], b_fract = b_fract)
#     # Save results
#     timesN1[Row] <- Sys.time() - start
#     file_name <- paste0(results_folder, "/ResultsN1Row", Row, ".RDS")
#     saveRDS(ssd_results, file = file_name)
# }


# ResultsN1: Results for determining n1
# ResultsN2: Results for determining n2


# Collect results in a big matrix ----------------------------------------------
## Fixed n1
all_results_N2 <- matrix(NA, ncol = 11, nrow_designN2 * b_fract)
row_result <- 1
for (row_desing in seq(nrow_designN2)) {
    readRDS(paste0("results/ResultsN1Row", row_desing, ".RDS"))
    row_result <- row_result
    for (b in seq(b_fract)) {
        median.BF01 <- median(ssd_results[[b]][[6]][, "BF.01"])
        median.BF10 <- median(ssd_results[[b]][[7]][, "BF.10"])
        mean.PMP0.H0 <- mean(ssd_results[[b]][[6]][, "PMP.0"])
        mean.PMP1.H0 <- mean(ssd_results[[b]][[6]][, "PMP.1"])
        mean.PMP0.H1 <- mean(ssd_results[[b]][[7]][, "PMP.0"])
        mean.PMP1.H1 <- mean(ssd_results[[b]][[7]][, "PMP.1"])
        n2 <- ssd_results[[b]]$n2
        eta.BF01 <- ssd_results[[b]]$Proportion.BF01
        eta.BF10 <- ssd_results[[b]]$Proportion.BF10
        n1 <- ssd_results[[b]]$n1
        all_results_N2[row_desing, ] <- c(b, median.BF01, median.BF10, mean.PMP0.H0,
                                          mean.PMP1.H0, mean.PMP0.H1, mean.PMP1.H1,
                                          n2, eta.BF01, eta.BF10, n1)
        row_result <- row_result + 1
    }
    
}
design_matrix_results <- design_matrixN2[rep(1:nrow(design_matrixN2), each = 3), ]
all_results_N2 <- as.data.frame(cbind(design_matrix_results, all_results_N2))
colnames(all_results_N2) <- c(names(design_matrixN2), "median.BF01", "median.BF10",
                              "mean.PMP0.H0", "mean.PMP1.H0", "mean.PMP0.H1",
                              "mean.PMP1.H1", "n2.final", "eta.BF01", "eta.BF10",
                              "n1.final")
saveRDS(all_results_N2, file = "results/all_results_N2.RDS")

## Fixed n2
all_results_N1 <- matrix(NA, ncol = 11, nrow_designN1 * b_fract)
row_result <- 1
for (row_desing in seq(nrow_designN1)) {
    readRDS(paste0("results/ResultsN1Row", row_desing, ".RDS"))
    row_result <- row_result
    for (b in seq(b_fract)) {
        median.BF01 <- median(ssd_results[[b]][[6]][, "BF.01"])
        median.BF10 <- median(ssd_results[[b]][[7]][, "BF.10"])
        mean.PMP0.H0 <- mean(ssd_results[[b]][[6]][, "PMP.0"])
        mean.PMP1.H0 <- mean(ssd_results[[b]][[6]][, "PMP.1"])
        mean.PMP0.H1 <- mean(ssd_results[[b]][[7]][, "PMP.0"])
        mean.PMP1.H1 <- mean(ssd_results[[b]][[7]][, "PMP.1"])
        n2 <- ssd_results[[b]]$n2
        eta.BF01 <- ssd_results[[b]]$Proportion.BF01
        eta.BF10 <- ssd_results[[b]]$Proportion.BF10
        n1 <- ssd_results[[b]]$n1
        all_results_N1[row_desing, ] <- c(b, median.BF01, median.BF10, mean.PMP0.H0,
                                          mean.PMP1.H0, mean.PMP0.H1, mean.PMP1.H1,
                                          n2, eta.BF01, eta.BF10, n1)
        row_result <- row_result + 1
    }
    
}
design_matrix_results <- design_matrixN1[rep(1:nrow(design_matrixN1), each = 3), ]
all_results_N1 <- as.data.frame(cbind(design_matrix_results, all_results_N1))
colnames(all_results_N1) <- c(names(design_matrixN1), "median.BF01", "median.BF10",
                              "mean.PMP0.H0","mean.PMP1.H0", "mean.PMP0.H1",
                              "mean.PMP1.H1", "n2.final", "eta.BF01", "eta.BF10",
                              "n1.final")
saveRDS(all_results_N1, file = "results/all_results_N1.RDS")


# all_results_N1: Results for determining n1
# all_results_N2: Results for determining n2

# Plots ------------------------------------------------------------------------