############################## Simulation for plots ###########################

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots

# Functions --------------------------------------------------------------------
source("SSD_clusters_null_plots.R")
#source("plots_SSD.R")

# Design matrix ----------------------------------------------------------------
n1 <- c(5, 10, 20, 40)
n2 <- 30
rho <- c(0.05, 0.1, 0.25) #Intraclass correlation
eff.size <- c(0.2, 0.5, 0.8)
bf.thresh <- 3
#fix <- c("n1", "n2") #Maybe not
design.matrix <- expand.grid(n1, n2, rho, eff.size, bf.thresh)
names(design.matrix) <- c("n1", "n2", "rho", "eff.size", "bf.thresh")

nrow.design <- nrow(design.matrix)
b <- 1
maxim <- 200
# Loop for every row -----------------------------------------------------------
for (Row in seq(nrow.design)) {
    # function
    ssd_results_null <- SSD_crt_null_plots(eff.size = design.matrix[Row, 4],
                                     n1 = design.matrix[Row, 1],
                                     n2 = design.matrix[Row, 2],
                                     n.datasets = 200,
                                     rho = design.matrix[Row, 3],
                                     BF.thresh = design.matrix[Row, 5],
                                     eta = 0.8, fixed = "n1", b.fract = b)
    # Save results
    save(ssd_results_null, file = paste("ResultPlotsRow", Row, ".Rdata", sep = ""))
    
}

# Collect results --------------------------------------------------------------

results_null_plots <- matrix(NA, ncol = 10, nrow = nrow.design * (86)) # To save the data for plots
r <- 0
for (i in seq(nrow.design)) {
    load(paste("ResultPlotsRow", i, ".Rdata", sep = ""))
    #browser()
    r <- r
    for (j in 2:87) {

        median.BF01 <- median(ssd_results_null[[j]][[6]][, "BF.01"]) # When null =TRUE
        median.BF10 <- median(ssd_results_null[[j]][[7]][, "BF.10"]) # When alternative = TRUE
        mean.PMP0.H0 <- mean(ssd_results_null[[j]][[6]][, "PMP.0"])## When null =TRUE
        mean.PMP1.H0 <- mean(ssd_results_null[[j]][[6]][, "PMP.1"])## When null =TRUE
        mean.PMP0.H1 <- mean(ssd_results_null[[j]][[7]][, "PMP.0"])## When alternative = TRUE
        mean.PMP1.H1 <- mean(ssd_results_null[[j]][[7]][, "PMP.1"])## When alternative = TRUE
        n1 <- ssd_results_null[[j]][["n1"]]
        n2 <- ssd_results_null[[j]][["n2"]] - 2
        eta.BF01 <- ssd_results_null[[j]][["Proportion.BF01"]]#
        eta.BF10 <- ssd_results_null[[j]][["Proportion.BF10"]]#
        r <-  r + 1
        results_null_plots[r, ] <- c(median.BF01, median.BF10, mean.PMP0.H0,
                                   mean.PMP1.H0, mean.PMP0.H1, mean.PMP1.H1,
                                   n1, n2, eta.BF01, eta.BF10)
        
    }
    
    
}
design_null <- design.matrix[ , -1]
design_null$ind <- seq(nrow.design)
design_null <- design_null[rep(1:nrow(design_null), times = 86), ]
design_null <- design_null[order(design_null$ind),]
results_all <- as.data.frame(cbind(design_null, results_null_plots))
head(results_null_plots)
head(results_all)
names(results_all) <- c(names(design_null), "med.BF01", "med.BF10",
                        "mean.PMP0.H0", "mean.PMP1.H0", "mean.PMP0.H1",
                        "mean.PMP1.H1", "n1", "n2.after", "eta.BF01", "eta.BF10")
save(results_all, file = "NullPlotsResults.Rdata")


# Plots --------------------------------
## Bayes factors --------------------
#rho and bf.threshold=same

ggplot(results_all, aes(med.BF01, color = as.factor(n1), fill = as.factor(n1))) +
    geom_histogram(alpha = 0.5, bins = 100) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Bayes Factor") + ylab("Frequency") #+ 
#scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

ggplot(results_all, aes(med.BF10, color = as.factor(n1), fill = as.factor(n1))) +
    geom_histogram(alpha = 0.5, bins = 100) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Bayes Factor") + ylab("Frequency") + 
    scale_y_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

##Final n2
ggplot(results_all, aes(y = med.BF01, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor")

ggplot(results_all, aes(y = med.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Bayes Factor H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor") +
    scale_y_log10(breaks = 10^(2:8), labels = trans_format("log10", math_format(10^.x)))# This is a problem
#This one doesn't work.
ggplot(results_all[results_all$eff.size == 0.2, ], aes(y = med.BF10, x = n2.after, color = as.factor(rho), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor") +
    scale_y_log10(breaks = 10^(2:8), labels = trans_format("log10", math_format(10^.x)))

ggplot(results_all[results_all$eff.size == 0.5, ], aes(y = med.BF10, x = n2.after, color = as.factor(rho), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor") +
    scale_y_log10(breaks = 10^(2:8), labels = trans_format("log10", math_format(10^.x)))

ggplot(results_all[results_all$eff.size == 0.8, ], aes(y = med.BF10, x = n2.after, color = as.factor(rho), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    labs(title = "Bayes Factor H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Bayes Factor") +
    scale_y_log10(breaks = 10^(2:4), labels = trans_format("log10", math_format(10^.x)))
## Eta ---------------------------------
ggplot(results_all, aes(y = eta.BF01, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Proportion BF H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
    geom_hline(yintercept = 0.8)

ggplot(results_all, aes(y = eta.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "Proportion BF H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
    geom_hline(yintercept = 0.8)

## PMP ------------------------------------------------------------------------ 
ggplot(results_all, aes(y = mean.PMP0.H0, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "PMP of H0 when H0 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Posterior Model Probabilities")

ggplot(results_all, aes(y = mean.PMP1.H1, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
    geom_point() + geom_line() +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
    labs(title = "PMP of H1 when H1 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
    xlab("Number of clusters") + ylab("Posterior Model Probabilities")

null_every <- results_all
