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
n2 <- 30
rho <- c(0.25, 0.1, 0.05) #Intraclass correlation
eff.size <- c(0.2, 0.5, 0.8)
bf.thresh <- c(1, 3, 5)
#fix <- c("n1", "n2") #Maybe not
design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh)
names(design.matrix) <- c("n1", "n2", "rho", "eff.size", "bf.thresh")

nrow.design <- nrow(design.matrix)
times_null <- matrix(NA, nrow = nrow.design, ncol = 1)
b <- 3

# Loop for every row -----------------------------------------------------------
for (Row in seq(nrow.design)) {
    start <- Sys.time()
    # function
    ssd_results_null <- SSD_crt_null(eff.size = design.matrix[Row, 4],
                                     n1 = design.matrix[Row, 1],
                                     n2 = design.matrix[Row, 2],
                                     n.datasets = 200,
                                     rho = design.matrix[Row, 3],
                                     BF.thresh = design.matrix[Row, 5],
                                     eta = 0.8, fixed = "n1", b.fract = b)
    # Save results
    times_null[Row] <- Sys.time() - start
    save(ssd_results_null, file = paste("ResultNullRow", Row, ".Rdata", sep = ""))
    
}
times_null <- as.data.frame(cbind(design.matrix, times_null))
save(times_null, file = "times_null.Rdata")
# Collect results --------------------------------------------------------------

results_null_all <- matrix(NA, ncol = 10, nrow = nrow.design * b)
r <- 1
for (i in seq(nrow.design)) {
    load(paste("ResultNullRow", i, ".Rdata", sep = ""))
    #load(file.path("ResultNullRow", i, ".Rdata", sep = ""))
    #browser()
    r <- r
    for (j in seq(b)) {
        median.BF01 <- median(ssd_results_null[[j]][[6]][, "BF.01"])
        median.BF10 <- median(ssd_results_null[[j]][[7]][, "BF.10"])
        mean.PMP0.H0 <- mean(ssd_results_null[[j]][[6]][, "PMP.0"])#
        mean.PMP1.H0 <- mean(ssd_results_null[[j]][[6]][, "PMP.1"])#
        mean.PMP0.H1 <- mean(ssd_results_null[[j]][[7]][, "PMP.0"])#
        mean.PMP1.H1 <- mean(ssd_results_null[[j]][[7]][, "PMP.1"])#
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
design_null <- design.matrix
design_null$ind <- seq(nrow.design)
design_null <- design_null[rep(1:nrow(design_null), times = 3), ]
design_null <- design_null[order(design_null$ind),]
b.frac <- rep(seq(b), nrow.design)
results_all <- as.data.frame(cbind(design_null, b.frac, results_null_all))
head(results_null_all)
head(results_all)
names(results_all) <- c(names(design_null), "b", "med.BF01", "med.BF10",
                             "mean.PMP0.H0", "mean.PMP1.H0", "mean.PMP0.H1",
                             "mean.PMP1.H1", "n2.final", "eta.BF01", "eta.BF10", "n1.final")
save(results_all, file = "AllNullResults.Rdata")



# Plots --------------------------------
## Bayes factors --------------------
#rho and bf.threshold=same
load("C:\\Users\\barra006\\OneDrive - Universiteit Utrecht\\Documents\\GitHub\\CRT-Sample-Size-Determination\\results\\output\\AllNullResults.Rdata")
load("AllNullResults.Rdata")
plots.SSD(1, data = results_all[which(results_all[ , "rho"] == 0.2 & results_all[ , "bf.thresh"] == 3), ], 
          y = med.BF01, grid_x = "n1.final", grid_y = "n2.final",
          title = "Bayes Factor H1 vs H2", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol",
          x_lab = "Bayes Factor", y_lab = "Frequency", color = eff.size)
plots.SSD(1, data = results_null_all, y = med.BF10, grid_x = n1, grid_y = n2,
          title = "Bayes Factor H1 vs H2", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol",
          x_lab = "Bayes Factor", y_lab = "Frequency", color = eff.size)
test <- results_all[which(results_all$bf.thresh == 3 & results_all$b == 1), ]
ggplot(test, aes(med.BF01, color = as.factor(n1), fill = as.factor(n1))) +
    geom_histogram(alpha = 0.5, bins = 100) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Bayes Factor") + ylab("Frequency") #+ 
    #scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

##Final n2
ggplot(test, aes(y = med.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor")

ggplot(test, aes(y = med.BF10, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor") +
    scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))# This is a problem
#Scale crazy!

## Eta ---------------------------------
ggplot(test, aes(y = eta.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Proportion BF H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Proportion of BF > BF threshold")

ggplot(test, aes(y = eta.BF10, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Proportion BF H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Proportion of BF > BF threshold")

## PMP ------------------------------------------------------------------------ 
ggplot(test, aes(y = mean.PMP0.H0, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "PMP of H0 when H0 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Posterior Model Probabilities")

ggplot(test, aes(y = mean.PMP1.H1, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "PMP of H1 when H1 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Posterior Model Probabilities")

## Sample size -------------------------
plots.SSD(2, data = results_all, x = n1, y = n2.final, grid_x = eff.size,
          grid_y = rho, title = "Numer of clusters in function of clusters' size",
          subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol", 
          x_lab = "Clusters' size",
          y_lab = "Number of clusters")


ggplot(test, aes(y = n2.final, x = n1, color = factor(n2), shape = factor(n2))) +
    geom_point() + geom_line() +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Numer of clusters in function of clusters' size", 
         subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters") 


### Final -------------------------------------------------------------------
ggplot(test[which(test$eff.size == 0.2),], 
       aes(y = n2.final, x = n1, color = factor(n2), shape = factor(n2))) +
    geom_point() + geom_line() + facet_grid(cols = vars(rho), labeller = label_both) +
    labs(title = "Number of clusters in function of clusters' size", 
         subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters")

ggplot(test[which(test$eff.size == 0.5),], 
       aes(y = n2.final, x = n1, color = factor(n2), shape = factor(n2))) +
    geom_point() + geom_line() + facet_grid(cols = vars(rho), labeller = label_both) +
    labs(title = "Number of clusters in function of clusters' size", 
         subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters")

ggplot(test[which(test$eff.size == 0.8),], 
       aes(y = n2.final, x = n1, color = factor(n2), shape = factor(n2))) +
    geom_point() + geom_line() + facet_grid(cols = vars(rho), labeller = label_both) +
    labs(title = "Number of clusters in function of clusters' size", 
         subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters")
# Maybe change the lines'width 
# Note that BF threshold and b fraction are constant!

########################### ONE CELL SIMULATION ################################
# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots
library(tictoc)

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
times_null <- matrix(NA, nrow = nrow.design, ncol = 1)
#times_null <- as.data.frame(cbind(design.matrix, times_null))
b <- 3
Row <- 1
start <- Sys.time()
ssd_results_null <- SSD_crt_null(eff.size = design.matrix[Row, 4],
                                 n1 = design.matrix[Row, 1],
                                 n2 = design.matrix[Row, 2],
                                 n.datasets = 100,
                                 rho = design.matrix[Row, 3],
                                 BF.thresh = design.matrix[Row, 5],
                                 eta = 0.8, fixed = "n1", b.fract = 1)
times_null[Row] <- Sys.time() - start
# Save results
save(ssd_results_null, file = paste("ResultNullRow", Row, ".Rdata", sep = ""))

#