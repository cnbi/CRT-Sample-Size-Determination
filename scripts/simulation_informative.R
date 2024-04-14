################################ SIMULATION INFORMATIVE HYPOTHESES ##################

# Libraries --------------------------------------------------------------------
library(dplyr) # For formatting the tables
library(ggplot2)
library(scales) # For scales in plots
library(latex2exp)


library(purrr) # For formatting the tables

# Functions --------------------------------------------------------------------
source("SSD_clusters_inform.R")

# Creation of folder for results
path <- "~/"
results_folder <- "resultsInform_FindN2"
dir.create(results_folder)

# Design matrix-----------------------------------------------------------------
# rho <- c()
rho <- c(0.025, 0.05, 0.1) #Change
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)

# Fixed n1
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(rho, eff_size, BF_threshold, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
colnames(design_matrixN2) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

#fixed n2
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")
nrow_designN1 <- nrow(design_matrixN1)
# ~N1: To determine n1
# ~N2: To determine n2

# Loop for every row -----------------------------------------------------------
## Find N2---------------------------------------------------------------------
for (Row in seq(nrow_designN2)) {
    # Start time
    start_time <- Sys.time()
    # Actual simulation
    ssd_results <- SSD_crt_inform(eff_size = design_matrixN2[Row, 2],
                                  n1 = design_matrixN2[Row, 4],
                                  ndatasets = 5000, rho = design_matrixN2[Row, 1],
                                  BF_thresh = design_matrixN2[Row, 3], eta = 0.8,
                                  fixed = as.character(design_matrixN2[Row, 5]),
                                  max = 1000, batch_size = 1000)
    # Save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0("ResultsN2Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN2Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    rm(ssd_results)
    gc()
}


## Find N1
for (Row in seq(nrow_designN1)) {
    # Start time
    start_time <- Sys.time()
    # Actual simulation
    ssd_results <- SSD_crt_null(eff_size = design_matrixN1[Row, 2],
                                n2 = design_matrixN1[Row, 4],
                                ndatasets = 5000, rho = design_matrixN1[Row, 1],
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                max = 500, batch_size = 1000)
    # Save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    rm(ssd_results)
    gc()
}
# ResultsN1: Results for determining n1
# ResultsN2: Results for determining n2

# Collect results in a big matrix ----------------------------------------------
## Find N2
all_results_N2 <- matrix(NA, ncol = 7, nrow = (nrow_designN2 ))
results_folder <- "resultsInform_FindN2"
# Big_result
for (row_design in seq(nrow_designN2)) {                                      
    stored_result <- readRDS(paste0(results_folder, "/ResultsN2Row", row_design, ".RDS")) # I am not sure if it's necessary to store the 
    
    median.BF12 <- median(stored_result[[4]][, "BF.12"])       
    median.BF21 <- median(stored_result[[4]][, "BF.21"])       
    mean.PMP1 <- mean(stored_result[[4]][, "PMP.1"])        
    mean.PMP2 <- mean(stored_result[[4]][, "PMP.2"])        
    eta.BF12 <- stored_result$Proportion.BF12
    n2 <- stored_result$n2
    n1 <- stored_result$n1
    all_results_N2[row_design, ] <- c(median.BF12, mean.PMP1, 
                                      median.BF21, mean.PMP2,
                                      eta.BF12, n2, n1)
}
all_results_N2 <- as.data.frame(cbind(design_matrixN2, all_results_N2))
colnames(all_results_N2) <- c(names(design_matrixN2), "median.BF12",
                              "mean.PMP1", "median.BF21", "mean.PMP2",
                              "eta.BF12", "n2.final", "n1.final")
saveRDS(all_results_N2, file = file.path(results_folder, "final_results_FindN2.RDS"))
# Running time in minutes
times_findN2 <- matrix(NA, nrow = nrow_designN2, ncol = 1)
for (row_result in seq(nrow_designN2)) {
    stored_result <- readRDS(paste0(results_folder, "/timeN2Row", row_result, ".RDS")) # I am not sure if it's necessary to store the 
    times_findN2[row_result, 1] <- stored_result
}
times_results_N2 <- as.data.frame(cbind(design_matrixN2, times_findN2))
colnames(times_results_N2) <- c(names(design_matrixN2), "total.time")             
saveRDS(times_results_N2, file = file.path(results_folder, "times_findN2.RDS"))


## Fixed n2
all_results_N1 <- matrix(NA, ncol = 7, nrow = nrow_designN1)
all_results_N1b <- matrix(NA, ncol = 7, nrow = nrow_designN1)
results_folder <- "resultsInform_FindN1-2"
for (row_design in seq(nrow_designN1)) {
    stored_result <- readRDS(paste0(results_folder, "/ResultsN1Row", row_design, ".RDS")) # I am not sure if it's necessary to store the 
    
    median.BF12 <- median(stored_result[[4]][, "BF.12"])       
    median.BF21 <- median(stored_result[[4]][, "BF.21"])       
    mean.PMP1 <- mean(stored_result[[4]][, "PMP.1"])        
    mean.PMP2 <- mean(stored_result[[4]][, "PMP.2"])        
    eta.BF12 <- stored_result$Proportion.BF12
    n2 <- stored_result$n2
    n1 <- stored_result$n1
    all_results_N1[row_design, ] <- c(median.BF12, mean.PMP1, 
                                      median.BF21, mean.PMP2,
                                      eta.BF12, n2, n1)
}
all_results_N1 <- as.data.frame(cbind(design_matrixN1, all_results_N1))
colnames(all_results_N1) <- c(names(design_matrixN1), "median.BF12",
                              "mean.PMP1", "median.BF21", "mean.PMP2",
                              "eta.BF12", "n2.final", "n1.final")
saveRDS(all_results_N1, file = file.path(results_folder, "final_results_FindN1.RDS"))
# Running time in minutes
times_findN1 <- matrix(NA, nrow = nrow_designN1, ncol = 1)
for (row_result in seq(nrow_designN1)) {
    stored_result <- readRDS(paste0(results_folder, "/timeN1Row", row_result, ".RDS")) # I am not sure if it's necessary to store the 
    times_findN1[row_result, 1] <- stored_result
}
times_results_N1 <- as.data.frame(cbind(design_matrixN1, times_findN1))
colnames(times_results_N1) <- c(names(design_matrixN1), "total.time")             
saveRDS(times_results_N1, file = file.path(results_folder, "final_times_findN1.RDS"))

# all_results_N1: Results for determining n1
# all_results_N2: Results for determining n2
