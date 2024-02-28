###################### FINAL SIMULATION FIRST PROJECT ##########################

# Libraries --------------------------------------------------------------------
library(dplyr) # For formatting the tables
library(ggplot2)
library(scales) # For scales in plots
library(latex2exp)


library(purrr) # For formatting the tables

# Functions --------------------------------------------------------------------
source("SSD_clusters_function.R")

# Creation of folder for results
path <- "~/"
results_folder <- "results"
dir.create(results_folder)

# Design matrix-----------------------------------------------------------------
# rho <- c()
rho <- c(0.025, 0.05, 0.1) #Change
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)
b_fract <- 3

# Fixed n1
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(rho, eff_size, BF_threshold, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
colnames(design_matrixN2) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

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

# ~N1: To determine n1
# ~N2: To determine n2

# Loop for every row -----------------------------------------------------------
## Find N2---------------------------------------------------------------------
for (Row in seq(nrow_designN2)) {
  # Start time
  start_time <- Sys.time()
  # Actual simulation
  ssd_results <- SSD_crt_null(eff_size = design_matrixN2[Row, 2],
                              n1 = design_matrixN2[Row, 4],
                              ndatasets = 5000, rho = design_matrixN2[Row, 1],
                              BF_thresh = design_matrixN2[Row, 3], eta = 0.8,
                              fixed = as.character(design_matrixN2[Row, 5]),
                              b_fract = b_fract, max = 1000, batch_size = 1000)
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
                              ndatasets = 10, rho = design_matrixN1[Row, 1],
                              BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                              fixed = as.character(design_matrixN1[Row, 5]),
                              b_fract = b_fract, max = 1000, batch_size = 100)
  # Save results
  end_time <- Sys.time()
  file_name <- file.path(results_folder, paste0( "serialResultsN1Row", Row, ".RDS"))
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
all_results_N2 <- matrix(NA, ncol = 11, nrow = (nrow_designN2 * b_fract))
row_result <- 1

# Big_result
for (row_design in seq(nrow_designN2)) {                                      
  stored_result <- readRDS(paste0(results_folder, "/ResultsN2Row", row_design, ".RDS")) # I am not sure if it's necessary to store the 
  row_result <- row_result
  for (b in seq(b_fract)) {
    median.BF01 <- median(stored_result[[b]][[6]][, "BF.01"])        # 6: data_H0
    median.BF10 <- median(stored_result[[b]][[7]][, "BF.10"])        # 7: data_H1
    mean.PMP0.H0 <- mean(stored_result[[b]][[6]][, "PMP.0"])        # 6: data_H0
    mean.PMP1.H0 <- mean(stored_result[[b]][[6]][, "PMP.1"])        # 6: data_H0
    mean.PMP0.H1 <- mean(stored_result[[b]][[7]][, "PMP.0"])        # 7: data_H1
    mean.PMP1.H1 <- mean(stored_result[[b]][[7]][, "PMP.1"])        # 7: data_H1
    n2 <- stored_result[[b]]$n2
    eta.BF01 <- stored_result[[b]]$Proportion.BF01
    eta.BF10 <- stored_result[[b]]$Proportion.BF10
    n1 <- stored_result[[b]]$n1
    all_results_N2[row_result, ] <- c(b, median.BF01, mean.PMP0.H0,
                                      mean.PMP1.H0, median.BF10, mean.PMP0.H1,
                                      mean.PMP1.H1, eta.BF01, eta.BF10, n2, n1)
    row_result <- row_result + 1
  }
}
design_matrix_results <- design_matrixN2[rep(1:nrow(design_matrixN2), each = 3), ]
all_results_N2 <- as.data.frame(cbind(design_matrix_results, all_results_N2))
colnames(all_results_N2) <- c(names(design_matrixN2), "b", "median.BF01",
                              "mean.PMP0.H0", "mean.PMP1.H0", "median.BF10",
                              "mean.PMP0.H1", "mean.PMP1.H1", "eta.BF01",
                              "eta.BF10", "n2.final", "n1.final")
saveRDS(all_results_N2, file = file.path(results_folder, "all_results_FindN2.RDS"))
# Running time in minutes
times_findN2 <- matrix(NA, nrow = length(1:40), ncol = 1)
for (row_result in seq(40)) {
  stored_result <- readRDS(paste0(results_folder, "/timeN2Row", row_result, ".RDS")) # I am not sure if it's necessary to store the 
  times_findN2[row_result, 1] <- stored_result
}
times_results_N2 <- as.data.frame(cbind(design_matrixN2[1:40, ], times_findN2))
colnames(times_results_N2) <- c(names(design_matrixN2), "total.time")             
saveRDS(times_results_N2, file = file.path(results_folder, "times_findN2_upto40.RDS"))



## Fixed n2
all_results_N1 <- matrix(NA, ncol = 11, nrow = (nrow_designN1 * b_fract))
row_result <- 1
for (row_desing in seq(nrow_designN1)) {
  stored_result <- readRDS(paste0(results_folder, "/ResultsN1Row", row_desing, ".RDS"))
  row_result <- row_result
  for (b in seq(b_fract)) {
    median.BF01 <- median(stored_result[[b]][[6]][, "BF.01"])
    median.BF10 <- median(stored_result[[b]][[7]][, "BF.10"])
    mean.PMP0.H0 <- mean(stored_result[[b]][[6]][, "PMP.0"])
    mean.PMP1.H0 <- mean(stored_result[[b]][[6]][, "PMP.1"])
    mean.PMP0.H1 <- mean(stored_result[[b]][[7]][, "PMP.0"])
    mean.PMP1.H1 <- mean(stored_result[[b]][[7]][, "PMP.1"])
    n2 <- stored_result[[b]]$n2
    eta.BF01 <- stored_result[[b]]$Proportion.BF01
    eta.BF10 <- stored_result[[b]]$Proportion.BF10
    n1 <- stored_result[[b]]$n1
    all_results_N1[row_desing, ] <- c(b, median.BF01, mean.PMP0.H0,
                                      mean.PMP1.H0, median.BF10, mean.PMP0.H1,
                                      mean.PMP1.H1, eta.BF01, eta.BF10, n2, n1)
    row_result <- row_result + 1
  }
}
design_matrix_results <- design_matrixN1[rep(1:nrow(design_matrixN1), each = 3), ]
all_results_N1 <- as.data.frame(cbind(design_matrix_results, all_results_N1))
colnames(all_results_N1) <- c(names(design_matrixN1), "b", "median.BF01",
                              "mean.PMP0.H0", "mean.PMP1.H0", "median.BF10",
                              "mean.PMP0.H1", "mean.PMP1.H1", "eta.BF01",
                              "eta.BF10", "n2.final", "n1.final")
saveRDS(all_results_N1, file = file.path(results_folder, "all_results_N1.RDS"))


# all_results_N1: Results for determining n1
# all_results_N2: Results for determining n2

# Plots ------------------------------------------------------------------------
## Bayes factors --------------------
#rho and bf.threshold=same

ggplot(all_results_N2, aes(median.BF01, color = as.factor(n1), fill = as.factor(n1))) +
  geom_histogram(alpha = 0.5, bins = 100) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = label_both) +
  labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Bayes Factor") + ylab("Frequency") #+ 
#scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

ggplot(all_results_N2, aes(median.BF10, color = as.factor(n1), fill = as.factor(n1))) +
  geom_histogram(alpha = 0.5, bins = 100) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = label_both) +
  labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Bayes Factor") + ylab("Frequency") + 
  scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

## Plot n2, Bayes factor ------------------------------------------------------------------------------------
## Change labels
rho_labs <- c("ICC: 0.025", "ICC: 0.05", "ICC: 0.1")
names(rho_labs) <- c("0.025", "0.05", "0.1")
eff_size_labs <- c(paste0("\u03B4 0.2"), paste0("\u03B4 0.5"), paste0("\u03B4 0.8"))
names(eff_size_labs) <- c("0.2", "0.5", "0.8")

ggplot(all_results_N2[all_results_N2$b == 1, ], aes(y = median.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = labeller(rho = rho_labs, eff_size = eff_size_labs)) +
  labs(title = "Bayes Factor H0 vs H1", color = "Cluster \nsizes", shape = "Cluster \nsizes") +
  xlab("Number of clusters") + ylab("Bayes Factor") + theme(legend.position = "bottom")

# Plot n2, Bayes factor comparing different b
ggplot(all_results_N2, aes(y = median.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line(aes(linetype = as.factor(b))) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = label_both) +
  labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Bayes Factor")

# Plot n2, Eta comparing different b
ggplot(all_results_N2[all_results$b == 1, ], aes(y = eta.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line(aes(linetype = as.factor(b))) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = label_both) +
  labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Eta")

# Save plot with high resolution
ggsave(filename = "n2BF", path = "results_SimulationFindN2", device='tiff', dpi=400)
ggsave(filename = "n2BF", path = "results_SimulationFindN2", width = 1070, height = 650, units = "px", device='png', dpi=400)


## Final n2: n2,------------------------------------------------------------------------------
ggplot(all_results_N2[all_results_N2$b == 1, ], aes(x = n1, y = n2.final, color = as.factor(BF_threshold), shape = as.factor(BF_threshold))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = labeller(rho = rho_labs, eff_size = eff_size_labs)) +
  labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nEffect sizes,and Intraclass Correlations",  
       color = "Bayes Factor \nThresholds", shape = "Bayes Factor \nThresholds") +
  xlab("Cluster sizes") + ylab("Number of clusters") + theme(legend.position = "bottom")

ggplot(all_results_N2, aes(x = n1, y = n2.final, color = as.factor(BF_threshold), shape = as.factor(b))) +
  geom_point() + geom_line(aes(linetype = as.factor(b))) +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff_size), labeller = labeller(rho = rho_labs, eff_size = eff_size_labs)) +
  labs(title = "Determining Number of Clusters in Function of Cluster Sizes, Bayes Factor Thresholds, \nEffect sizes,and Intraclass Correlations",  
       color = "Bayes Factor \nThresholds", shape = "Fraction b", linetype = "Fraction b") +
  xlab("Cluster sizes") + ylab("Number of clusters") + theme(legend.position = "bottom")

# TABLES -------------------------------------------------------------------------------------------------
## Option 1
library(xtable)
table_results <- c("rho", "BF_threshold", "n1", "n2.final", "eta.BF10", "eta.BF01")
round_col <- c("BF_threshold", "n1", "n2.final")
results_colnames <- c("ICC", "$BF_{thresh}$", "N1", "N2", "$P(BF_{10}>BF_{thresh})$", "$P(BF_{01}>BF_{thresh})$")
effect.size0.2 <- all_results_FindN2[which(all_results_FindN2$b == 1 & all_results_FindN2$eff_size == 0.2),table_results]
effect.size0.2 <- effect.size0.2[order(effect.size0.2$rho, effect.size0.2$BF_threshold, effect.size0.2$n1), ]
effect.size0.2[, round_col] <- round(effect.size0.2[, round_col], digits = 0)
colnames(effect.size0.2) <- results_colnames
decimals <- c(3, 0, 0, 0, 3, 3)
for (i in 1:length(table_results)) {
  effect.size0.2[, i] <- formatC(effect.size0.2[, i], format = "f", digits = decimals[i])
}
print(xtable(effect.size0.2), include.rownames = FALSE)

effect.size0.5 <- all_results_FindN2[which(all_results_FindN2$b == 1 & all_results_FindN2$eff_size == 0.5),table_results]
effect.size0.8 <- all_results_FindN2[which(all_results_FindN2$b == 1 & all_results_FindN2$eff_size == 0.8),table_results]

# Option 2
#rownames(try_tr) <- NULL #To reset the index
try_tr <- effect.size0.2[rep(1:nrow(effect.size0.2), each = 2), ]
try_tr[1:nrow(try_tr) %% 2 == 0, ] <- NA
try_tr <- try_tr[, 1:ncol(effect.size0.2) - 1]
colnames(try_tr) <- c("ICC", "$BF_{thresh}$", "N1", "N2", "$P(BF_>BF_{thresh})$")
try_tr[1:nrow(try_tr) %% 2 == 0, "$P(BF_>BF_{thresh})$"] <- effect.size0.2$`$P(BF_{01}>BF_{thresh})$`
try_tr$scenario <- c("H1", "H0")
nrow_result <- nrow(try_tr)/3
results_new_format <- try_tr[1:nrow_result, 2:ncol(try_tr)]
for (index in 1:2) {
  new_col <- try_tr[((nrow_result*index) + 1):((nrow_result*index) + nrow_result), c("N2", "$P(BF_>BF_{thresh})$")]
  browser()
  results_new_format <- cbind(results_new_format, new_col)
}
results_new_format <- results_new_format[, c(1, 2, 5, 3, 4, 6, 7, 8, 9)]
print(xtable(results_new_format), include.rownames = FALSE)


# Find n1
table_results <- c("rho", "BF_threshold", "n2", "n1.final", "eta.BF10", "eta.BF01")
effect.size0.2N1 <- all_results_findN1[which(all_results_findN1$b == 1 & all_results_findN1$eff_size == 0.2),table_results]
