########################## TEST NULL HYPOTHESIS ################################

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots

# Functions --------------------------------------------------------------------
source("SSD_clusters_function.R")
source("plots_SSD.R")

# Design matrix ----------------------------------------------------------------
n1 <- c(5, 10, 20, 40)
n2 <- c(30, 60 ,90)
rho <- c(0.25, 0.05, 0.1) #Intraclass correlation
eff.size <- c(0.2, 0.5, 0.8)
bf.thresh <- c(1, 3, 5)
#fix <- c("n1", "n2") #Maybe not
design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh)
names(design.matrix) <- c("n1", "n2", "rho", "eff.size", "bf.thresh")

nrow.design <- nrow(design.matrix)
times_null <- matrix(NA, nrow = nrow.design, ncol = 13)

b <- 3

# Loop for every row -----------------------------------------------------------
for (Row in seq(nrow.design)) {
    # function
    times_null[Row, ] <- bench::mark( ssd_results_null <- SSD_crt_null(eff.size = design.matrix[Row, 4],
                                     n1 = design.matrix[Row, 1],
                                     n2 = design.matrix[Row, 2],
                                     n.datasets = 1000,
                                     rho = design.matrix[Row, 3],
                                     BF.thresh = design.matrix[Row, 5],
                                     eta = 0.8, fixed = "n1", b.fract = b) )
    # Save results
    save(ssd_results_null, file = paste("ResultNullRow", Row, ".Rdata", sep = ""))
    
}
save(times_null, file = "times_null.Rdata")
# Collect results --------------------------------------------------------------

results_null_all <- matrix(NA, ncol = 6, nrow = nrow.design * b)
r <- 1
for (i in seq(nrow.design)) {
    load(file.path("ResultsRow", i, ".Rdata", sep = ""))
    r <- r
    for (j in seq(b)) {
        median.BF01 <- median(ssd_results_null[[j]][[6]][, "BF.01"])
        median.BF10 <- median(ssd_results_null[[j]][[7]][, "BF.10"])
        mean.PMP0.H0 <- mean(ssd_results_null[[j]][[6]][, "BF.01"])#
        mean.PMP1.H0 <- mean(ssd_results_null[[j]][[6]][, "BF.01"])#
        mean.PMP0.H1 <- mean(ssd_results_null[[j]][[7]][, "BF.01"])#
        mean.PMP1.H1 <- mean(ssd_results_null[[j]][[7]][, "BF.01"])#
        n2 <- ssd_results_null[[j]]$n2
        eta.BF01 <- ssd_results_null[[j]]$Proportion.BF01#
        eta.BF10 <- ssd_results_null[[j]]$Proportion.BF10#
        n1 <- ssd_results_null[[j]]$n1
        results_null_all[r, ] <- c(median.BF01, median.BF10, mean.PMP0.H0,
                                   mean.PMP1.H0, mean.PMP0.H1, mean.PMP1.H1,
                                   n2, eta.BF01, eta.BF10, n1)
        r = r + 1
    }
    
    
}
design_null <- rep(design.matrix, each = b)
b.frac <- rep(seq(b), nrow.design)
results_all <- as.data.frame(cbind(design_null, b.frac, results_null_all))
head(results_null_all)
names(results_null_all) <- c(names(design.matrix), b.frac, "med.BF01", "med.BF10",
                             "mean.PMP0.H0", "mean.PMP1.H0", "mean.PMP0.H1",
                             "mean.PMP1.H1", "n2", "eta.BF01", "eta.BF10", "n1")
save(results_null_all, file = "AllNullResults.Rdata")

# Plots --------------------------------
## Bayes factors --------------------

plots.SSD(1, data = results_null_all, y = med.BF01, grid_x = n1, grid_y = n2,
          title = "Bayes Factor H1 vs H2", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol",
          x_lab = "Bayes Factor", y_lab = "Frequency", color = eff.size)
plots.SSD(1, data = results_null_all, y = med.BF10, grid_x = n1, grid_y = n2,
          title = "Bayes Factor H1 vs H2", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol",
          x_lab = "Bayes Factor", y_lab = "Frequency", color = eff.size)
## Sample size -------------------------
plots.SSD(2, data = results_null_all, x = n1, y = mean.n2, grid_x = eff.size,
          grid_y = rho, title = "Numer of clusters in function of clusters' size",
          subtitle = "H1:Dintervention>Dcontrol", x_lab = "Clusters' size",
          y_lab = "Number of clusters")