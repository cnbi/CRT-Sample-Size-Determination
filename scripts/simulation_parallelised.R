


#Libraries
install.packages("doParallel")
install.packages("foreach")
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)
library(ps)

source("SSD_clusters_function.R")
# Creation of folder for results
path <- "~/"
results_folder <- "results"

# Design matrix-----------------------------------------------------------------
rho <- c(0.025, 0.05, 0.1) #Change
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)
b_fract <- 3
#fixed n2
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")
# Rows that I already ran find N1 
rho <- c(0.05, 0.1)
eff_size <- c(0.2)
BF_threshold <- c(1, 3, 5)
b_fract <- 3
n2 <- c(30, 60, 90)
completed_findN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
colnames(completed_findN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")

# Filter Find N1
completed_findN1 <- as.data.frame(completed_findN1)
design_matrixN1 <- as.data.frame(design_matrixN1)
design_matrixN1 <- anti_join(design_matrixN1, completed_findN1)
nrow_designN1 <- nrow(design_matrixN1)

# Simulation parallelised---------------------------------------------------------
# Detect
ncluster <- detectCores()/2
ncluster <- 2
# Create clusters and register them
cl <- makeCluster(ncluster)
registerDoParallel(cl)
# Export libraries, functions and variables
clusterExport(cl, c("SSD_crt_null", "results_folder", "design_matrixN1", "ps", "b_fract"))
# Divide rows in clusters
identify_clusters <- cut((15:nrow_designN1), breaks = ncluster, labels = FALSE)
rows_divided <- split(15:nrow_designN1, 1:ncluster)
rows_divided <- lapply(15:nrow_designN1, function(x) list(x))
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
                                ndatasets = 100, rho = design_matrixN1[Row, 1], #change
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                b_fract = b_fract, max = 1000, batch_size = 1000)
    # Save results
    end_time <- Sys.time()
    eff_size = design_matrixN1[Row, 2]
    is.vector(Row)
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
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
ncluster <- detectCores()/2
ncluster <- 2
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
                                ndatasets = 10, rho = design_matrixN1[Row, 1], #change
                                BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                fixed = as.character(design_matrixN1[Row, 5]),
                                b_fract = b_fract, max = 1000, batch_size = 1000)
    # Save results
    end_time <- Sys.time()
    eff_size = design_matrixN1[Row, 2]
    is.vector(Row)
    file_name <- file.path(results_folder, paste0("ResultsN1Row", Row, ".RDS"))
    saveRDS(ssd_results, file = file_name)
    # Save running time
    running_time <- as.numeric(difftime(end_time,start_time, units = "mins"))
    time_name <- file.path(results_folder, paste0("timeN1Row", Row, ".RDS"))
    saveRDS(running_time, file = time_name)
    rm(ssd_results)
    gc()
    }
}

foreach(Rows = rows_divided) %dopar% {
    run_simulation(Rows)
}

stopCluster(cl)
stopImplicitCluster()
