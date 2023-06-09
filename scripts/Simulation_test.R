###########################TEST INFORMATIVE HYPOTHESES ############################

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots


# Functions needed -------------------------------------------------------------
source("SSD_clusters_inform.R")
source("plots_SSD.R")

# Design matrix ----------------------------------------------------------------
n1 <- c(5, 10, 20, 40)
n2 <- c(30, 60 ,90)
rho <- c(0.25, 0.05, 0.1) #Intraclass correlation
eff.size <- c(0.2, 0.5, 0.8)
bf.thresh <- c(1, 3, 5)
#fix <- c("n1", "n2")
#design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh, fix)
design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh)
#names(design.matrix) <- c("n1", "n2", "rho", "eff.size", "bf.thresh", "fix")
names(design.matrix) <- c("n1", "n2", "rho", "eff.size", "bf.thresh")

nrow.design <- nrow(design.matrix)
times <- matrix(NA, nrow = nrow.design, ncol = 13) 

# Loop for every row -----------------------------------------------------------
for (Row in seq(nrow.design)) {
    start <- Sys.time()
    # function
    ssd_results <- SSD_crt_inform(eff.size = design.matrix[Row, 4], 
                                                  n1 = design.matrix[Row, 1],
                                                  n2 = design.matrix[Row, 2],
                                                  n.datasets = 1000,
                                                  rho = design.matrix[Row, 3],
                                                  BF.thresh = design.matrix[Row, 5],
                                                  eta = 0.8, 
                                                  fixed = "n1")

    # Save results
    times[Row] <- Sys.time() - start
    save(ssd_results, file = paste("ResultsRow", Row, ".Rdata", sep = ""))
}
times <- as.data.frame(cbind(design.matrix, times))
save(times, file = "times_inform.Rdata")

# Collect results --------------------------------------------------------------
results_all <- matrix(NA, ncol = 6, nrow = nrow.design)
for (i in seq(nrow.design)) {
    load(file.path("ResultsRow", i, ".Rdata", sep = ""))
    median.BF12 <- median(ssd_results[[4]][, "BF.12"])
    median.BF21 <- median(ssd_results[[4]][, "BF.21"]) # I think I can eliminate this one
    mean.PMP1 <- mean(ssd_results[[4]][, "PMP.1"])
    mean.PMP2 <- mean(ssd_results[[4]][, "PMP.2"])
    n2 <- ssd_results$n2
    eta <- ssd_results$eta
    n1 <- ssd_results$n1
    results_all[i, ] <- c(median.BF21, median.BF12, mean.PMP1, mean.PMP2, n2,
                          eta, n1)
}
head(results_all)
#names(design_results) <- names(design.matrix)
results_all <- as.data.frame(cbind(design.matrix, results_all))
names(results_all) <- c(names(design.matrix), "med.BF12", "med.BF21", "mean.PMP1",
                        "mean.PMP2", "n2", "eta", "n1")
save(results_all, file = "AllResults.Rdata")

# Plots ------------------------------------------------------------------------------
## Bayes factors -------------------------------------
hist(log10(results_all$median.BF12), breaks = 100)
ggplot(results_all, aes(median.BF12, color = eff.size, fill = eff.size)) +
    geom_histogram(alpha = 0.5, bins = 50) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(n1), cols = vars(n2), labeller = label_both) +
    labs(title = "Bayes Factor H1 vs H2", subtitle = "H1:Dintervention>Dcontrol") +
    xlab("Bayes Factor") + ylab("Frequency") + 
    scale_x_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

plots.SSD(plot = 1, data = results_all, x = median.BF12, grid_x = n1, grid_y = n2,
          title = "Bayes Factor H1 vs H2", subtitle = "H1:Dintervention>Dcontrol",
          x_lab = "Bayes Factor", y_lab = "Frequency", color = eff.size)
## sample size ---------------------------------------
ggplot(results_grouped, aes(y = mean.n2, x = n1)) + geom_line() +
    facet_grid(rows = vars(eff.size), cols = vars(rho), labeller = label_both) +
    labs(title = "Numer of clusters in function of clusters' size", 
         subtitle = "H1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters") 

plots.SSD(plot = 2, data = results_all, x = n1, y = mean.n2, grid_x = eff.size,
          grid_y = rho, title = "Numer of clusters in function of clusters' size",
          subtitle = "H1:Dintervention>Dcontrol", x_lab = "Clusters' size",
          y_lab = "Number of clusters")