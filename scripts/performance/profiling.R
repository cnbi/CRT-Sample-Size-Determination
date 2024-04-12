############################## PROFILING ###############################

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots
library(profvis)
source("SSD_clusters_function.R")

# Design matrix
rho <- c(0.05, 0.1) #Change
eff_size <- c(0.2)
BF_threshold <- c(1, 3, 5)
b_fract <- 3
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
nrow_designN1 <- nrow(design_matrixN1)
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n2", "fixed")

profvis(
    # run simulation
    for (Row in seq(nrow_designN1)) {
        # Start time
        startt <- Sys.time()
        # Actual simulation
        ssd_results <- SSD_crt_null(eff_size = design_matrixN1[Row, 2],
                                    n1 = design_matrixN1[Row, 4],
                                    ndatasets = 10, rho = design_matrixN1[Row, 1],
                                    BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                    fixed = as.character(design_matrixN1[Row, 5]),
                                    b_fract = b_fract)
        # Save results
        endd <- Sys.time()
        file_name <- file.path(paste0("profiling_row", Row, ".RDS"))
        saveRDS(ssd_results, file = file_name)
        diff.time <- as.numeric(difftime(endd, startt, units = "mins"))
        time_name <- file.path(paste0("profiling_time_row", Row, ".RDS"))
        saveRDS(diff.time, file = time_name)
    }
)