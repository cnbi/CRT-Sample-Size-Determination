################################# PARALLELISED VERSION ########################

#Libraries
if (!require(parallel)) {install.packages("parallel")}
if (!require(foreach)) {install.packages("foreach")}
if (!require(doParallel)) {install.packages("doParallel")}
if (!require(dplyr)) {install.packages("dplyr")}

# Function
source("SSD_clusters_inform.R")

# Creation of folder for results
path <- "~/"
results_folder <- "resultsInform_FindN2"
dir.create(results_folder)

results_folder <- "resultsInform_FindN1"
dir.create(results_folder)

# Design matrix-----------------------------------------------------------------
rho <- c(0.025, 0.05, 0.1)
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)

# Fixed n1: Finding n2
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(rho, eff_size, BF_threshold, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
colnames(design_matrixN2) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

# Fixed n2: Finding n1
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")
nrow_designN1 <- nrow(design_matrixN1)

# ~N1: To determine n1
# ~N2: To determine n2

# Foreach ---------------------------------------------------------------------
# Detect clusters
# ncluster <- detectCores()/2
ncluster <- 3

# Create clusters and register them
cl <- makeCluster(ncluster)
registerDoParallel(cl)
rows_divided <- split(1:nrow_designN1, 1:ncluster)
clusterExport(cl, c("SSD_crt_inform", "results_folder", "design_matrixN1"))

# n2 fixed: Find n1
results_folder <- "resultsInform_FindN1"
run_simulation <- function(Rows) {
  for (Row in Rows) {
    # Start time
    start_time <- Sys.time()
    
    # Actual simulation
    ssd_results <- SSD_crt_inform(eff_size = design_matrixN1[Row, 2],
                                n2 = design_matrixN1[Row, 4],
                                ndatasets = 5000, rho = design_matrixN1[Row, 1],
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                max = 1000, batch_size = 1000)
    
    # Save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    
    # Save running time
    running_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    # Clean
    rm(ssd_results)
    gc()
  }
}

#n1 fixed: Find n2
results_folder <- "resultsInform_FindN2"
run_simulation <- function(Rows) {
  for (Row in Rows) {
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
    running_time <- as.numeric(difftime(end_time, start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN2Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    # Clean
    rm(ssd_results)
    gc()
  }
}


foreach(Rows = rows_divided) %dopar% {
  run_simulation(Rows)
}

# Stop clusters
stopCluster(cl)
stopImplicitCluster()
