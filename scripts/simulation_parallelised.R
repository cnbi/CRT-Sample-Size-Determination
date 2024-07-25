############## Parallelisation of simulation for null hypothesis ###############


# Libraries
if (!require(parallel)) {install.packages("parallel")}
if (!require(foreach)) {install.packages("foreach")}
if (!require(doParallel)) {install.packages("doParallel")}
if (!require(dplyr)) {install.packages("dplyr")}
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(ps)

# Function
source("SSD_clusters_function.R")

# Creation of folder for results
path <- "~/"
results_folder <- "resultsNull"

# Design matrix-----------------------------------------------------------------
rho <- c(0.025, 0.05, 0.1)
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)
b_fract <- 3

# Fixed n2
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")
nrow_designN1 <- nrow(design_matrixN1)

# Simulation parallelised-------------------------------------------------------
# Detect
ncluster <- detectCores() / 2
#ncluster <- 2

# Create clusters and register them
cl <- makeCluster(ncluster)
registerDoParallel(cl)

# Export libraries, functions and variables
clusterExport(cl, c("SSD_crt_null", "results_folder", "design_matrixN1", "ps", "b_fract"))

# Divide rows in clusters
identify_clusters <- cut(seq(nrow_designN1), breaks = ncluster, labels = FALSE)
rows_divided <- split(seq(nrow_designN1), 1:ncluster)
rows_divided <- lapply(seq(nrow_designN1), function(x) list(x))

# Creation of simulation function
run_simulation <- function(Row) {
    # Start time
    start_time <- Sys.time()

    # Status
    cat("Running row", Row, "\n")
    resources <- ps()
    print(resources)

    # Actual simulation
    ssd_results <- SSD_crt_null(eff_size = design_matrixN1[Row, 2],
                                n2 = design_matrixN1[Row, 4],
                                ndatasets = 1000, rho = design_matrixN1[Row, 1],
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                b_fract = b_fract, max = 1000, batch_size = 1000)

    # Save results
    end_time <- Sys.time()
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    
    # Clean
    rm(ssd_results)
    gc()
}

# Apply simulation
clusterApply(cl, rows_divided, run_simulation)

# Stop parallelisation
stopCluster(cl)
stopImplicitCluster()

# Foreach ---------------------------------------------------------------------
# Detect
ncluster <- detectCores() / 2
# ncluster <- 2

# Create clusters and register them
cl <- makeCluster(ncluster)
registerDoParallel(cl)
rows_divided <- split(15:nrow_designN1, 1:ncluster)
clusterExport(cl, c("SSD_crt_null", "results_folder", "design_matrixN1", "ps", "b_fract"))

run_simulation <- function(Rows) {
    for (Row in Rows) {
    # Start time
    start_time <- Sys.time()

    # Status
    cat("Running row", Row, "\n")
    resources <- ps()
    print(resources)

    # Actual simulation
    ssd_results <- SSD_crt_null(eff_size = design_matrixN1[Row, 2],
                                n2 = design_matrixN1[Row, 4],
                                ndatasets = 1000, rho = design_matrixN1[Row, 1],
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                b_fract = b_fract, max = 1000, batch_size = 1000)

    # Save results
    end_time <- Sys.time()
    is.vector(Row)
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)

    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
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
