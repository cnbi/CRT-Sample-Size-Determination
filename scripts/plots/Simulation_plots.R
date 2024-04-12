############################## Simulation for plots ###########################

# Libraries --------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots

# Functions --------------------------------------------------------------------
source("SSD_clusters_function_pl.R")
#source("plots_SSD.R")

# Creation of folder for results
path <- "~/"
results_folder <- "results"
dir.create(results_folder)

# Design matrix ----------------------------------------------------------------
rho <- c(0.025, 0.05, 0.1)
eff_size <- c(0.2, 0.5, 0.8)
BF_threshold <- c(1, 3, 5)
b_fract <- 3
# fixed <- c("n1", "n2")
# Fixed n1
n1 <- c(5, 10, 40)
design_matrixN2 <- expand.grid(rho, eff_size, BF_threshold, n1, fixed <- "n1")
nrow_designN2 <- nrow(design_matrixN2)
# timesN2 <- matrix(NA, nrow = nrow_designN2, ncol = 1)
timesN2_1 <- matrix(NA, nrow = length(41:81), ncol = 1)
colnames(design_matrixN2) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

#fixed n2
n2 <- c(30, 60, 90)
design_matrixN1 <- expand.grid(rho, eff_size, BF_threshold, n2, fixed <- "n2")
nrow_designN1 <- nrow(design_matrixN1)
# timesN1 <- matrix(NA, nrow = nrow_designN1, ncol = 1)
colnames(design_matrixN1) <- c("rho", "eff_size", "BF_threshold", "n1", "fixed")

# ~N1: To determine n1
# ~N2: To determine n2
# Loop for every row -----------------------------------------------------------
# Finding n2 (i.e. fixed n1)
for (Row in 1:40) {
  # Start time
  startt <- Sys.time()
  # Actual simulation
  ssd_results <- SSD_crt_null_plots(eff_size = design_matrixN2[Row, 2],
                                    n1 = design_matrixN2[Row, 4],
                                    ndatasets = 500, rho = design_matrixN2[Row, 1],
                                    BF_thresh = design_matrixN2[Row, 3], eta = 0.8,
                                    fixed = as.character(design_matrixN2[Row, 5]), b_fract = b_fract)
  # Save results
  endd <- Sys.time()
  
  file_name <- file.path(paste0(results_folder, "/ResultsN2Row", Row, ".RDS"))
  saveRDS(ssd_results, file = file_name)
  #Save running time
  timesN2_1[Row] <- as.numeric(difftime(endd, startt, units = "mins"))
  
}

saveRDS(timesN2_1, file = file.path("results_plots\timesN2_.RDS"))
# Finding n2 (i.e. fixed n1)
for (Row in 1:40) {
  # Start time
  # Start time
  start <- Sys.time()
  # Actual simulation
  ssd_results <- SSD_crt_null_plots(eff_size = design_matrixN1[Row, 2],
                                    n2 = design_matrixN1[Row, 4],
                                    ndatasets = 5000, rho = design_matrixN1[Row, 1],
                                    BF_thresh = design_matrixN1[Row, 3], eta = 0.8,
                                    fixed = design_matrixN1[Row, 5], b_fract = b_fract)
  # Save results
  # timesN1_row <- Sys.time() - start
  file_name <- file.path(paste0(results_folder, "/ResultsN1Row", Row, ".RDS"))
  saveRDS(ssd_results, file = file_name)
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
results16 <- load("~/Results-19062023/NullPlotsResults.Rdata")
results25 <- load("~/Results-25062023/AllNullResults.Rdata")

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
rho_labs <- c("ICC: 0.025", "ICC: 0.05", "ICC: 0.1")
names(rho_labs) <- c("0.025", "0.05", "0.1")
eff_size_labs <- c(paste0("\u03B4 0.2"), paste0("\u03B4 0.5"), paste0("\u03B4 0.8"))
names(eff_size_labs) <- c("0.2", "0.5", "0.8")
ggplot(results_all[results_all$rho != 0.25,], aes(y = eta.BF01, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), labeller = labeller(rho = rho_labs)) +
  xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
  labs(color = "Cluster \nsizes", shape = "Cluster \nsizes") +
  geom_hline(yintercept = 0.8, colour = "red", linetype = "dashed") + theme(legend.position = "bottom")
ggsave(filename = "n2Eta01.png", path = "C:/Users/barra006/OneDrive - Universiteit Utrecht/Documents/GitHub/CRT-Sample-Size-Determination/scripts", width = 1070, height = 650, units = "px", device='png', dpi=400)

ggplot(results_all[results_all$rho != 0.25,], aes(y = eta.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = labeller(rho = rho_labs, eff.size = eff_size_labs)) +
  xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
  labs(color = "Cluster \nsizes", shape = "Cluster \nsizes") +
  geom_hline(yintercept = 0.8, colour = "red", linetype = "dashed") + theme(legend.position = "bottom")
ggsave(filename = "n2Eta10.png", path = "C:/Users/barra006/OneDrive - Universiteit Utrecht/Documents/GitHub/CRT-Sample-Size-Determination/scripts", width = 1070, height = 650, units = "px", device='png', dpi=400)

ggplot(results_all, aes(y = eta.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), labeller = label_both) +
  labs(title = "Proportion BF H0 vs H1", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
  geom_hline(yintercept = 0.8, colour = "red", linetype = "dashed")

rho.labs <- c("ICC: 0.05", "ICC: 0.1", "ICC: 0.25")
names(rho.labs) <- c("0.05", "0.1", "0.25")
tiff("propBF01.tif", width = 3100, height = 1800, units = "px", res = 300)
ggplot(results_all[results_all$eff.size == 0.2,], aes(y = eta.BF01, x = n2.final, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(cols = vars(rho), labeller = labeller(rho = rho.labs)) +
  labs(title = "Proportion of Bayes Factors Exceeding the Threshold under the Null Hypothesis", subtitle = "Effect size: 0.2",
       color = "Cluster sizes", shape = "Cluster sizes") +
  xlab("Number of clusters") + ylab("Proportion of BF01") +
  geom_hline(yintercept = 0.8) 
dev.off()



ggplot(results_all, aes(y = eta.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
  labs(title = "Proportion BF H1 vs H0", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Proportion of BF > BF threshold") +
  geom_hline(yintercept = 0.8)

tiff("propBF10.tif", width = 3100, height = 1800, units = "px", res = 300)
ggplot(results_all[which(results_all$eff.size == 0.2 & results_all$bf.thresh == 3),], aes(y = eta.BF10, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(cols = vars(rho), labeller = labeller(rho = rho.labs)) +
  labs(title = "Proportion of Bayes Factors Exceeding the Threshold under the Alternative Hypothesis", subtitle = "Effect size: 0.2",
       color = "Cluster sizes", shape = "Cluster sizes") +
  xlab("Number of clusters") + ylab("Proportion of BF10") +
  geom_hline(yintercept = 0.8)# + labs(color = "Cluster sizes", shape = "Cluster sizes")
dev.off()
## PMP ------------------------------------------------------------------------ 
ggplot(results_all, aes(y = mean.PMP0.H0, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
  labs(title = "PMP of H0 when H0 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dihttp://127.0.0.1:44319/graphics/plot_zoom_png?width=1200&height=778ntervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Posterior Model Probabilities")

ggplot(results_all, aes(y = mean.PMP1.H1, x = n2.after, color = as.factor(n1), shape = as.factor(n1))) +
  geom_point() + geom_line() +
  scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  facet_grid(rows = vars(rho), cols = vars(eff.size), labeller = label_both) +
  labs(title = "PMP of H1 when H1 = TRUE", subtitle = "H0:Dintervention=Dcontrol \nH1:Dintervention>Dcontrol") +
  xlab("Number of clusters") + ylab("Posterior Model Probabilities")

null_every <- results_all


# Sampling plot --------------------------------------
library(ggplot2)
sampling5000 <- SSD_crt_null_plots(eff_size = 0.5, n1 = 10, n2 = 30, ndatasets = 5000, rho = 0.03, 
                   BF_thresh = 3, eta = 0.8, fixed = "n2", b_fract = 3, increasing = FALSE,
                   max = 100, batch_size = 500)

#b=1
dataH0 <- as.data.frame(sampling5000[[1]][6])
names(dataH0)
ggplot(dataH0, aes(x = data_H0.BF.01)) + geom_density(alpha = 0.3) + xlim(0, 50) +
  geom_vline(aes(xintercept = 3),linetype = "dashed", linewidth = 1, colour = "red")
polygon(c(dataH0[dataH0$data_H0.BF.01 > 3], max(dataH0$data_H0.BF.01), 3), c(density[dataH0$data_H0.BF.01 >3], 0, 0, col = "grey"))

densities <- density(dataH0$data_H0.BF.01)
d <- data.frame(x = densities$x, y = densities$y) # I think this is the same
dd <- with(densities,data.frame(x,y))  # I think this is the same
qplot(x,y,data=dd,geom="line") +
  geom_ribbon(data=subset(dd,x>3 & x<50),aes(ymax=y),ymin=0,
              fill="red",colour=NA,alpha=0.5)
x1 <- min(which(densities$x >= 3))  
x2 <- max(which(densities$x <  50))
with(densities, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col="gray"))
