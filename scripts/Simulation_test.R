###########################TEST INFORMATIVE HYPOTHESES ############################

# Libraries
library(dplyr)
library(ggplot2)
library(scales) # For scales in plots

# Design matrix
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

# Loop for every row
for (Row in seq(nrow.design)) {
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
    save(ssd_results, file = paste("ResultsRow", Row, ".Rdata", sep = ""))

}

# Collect results
nrow.results <- nrow(ssd_results[[4]])
results_all <- matrix(NA, ncol = ncol(ssd_results[[4]]), nrow = nrow.results * nrow.design)
for (i in seq(nrow.design)) {
    Row <- i
    load(file.path("ResultsRow", Row, ".Rdata", sep = ""))
    results_all[nrow.results * (i - 1) + 1:(i * nrow.results), ]
}
head(results_all)
design_results <- design.matrix[rep(1:nrow(design.matrix),each = nrow.results),]
#names(design_results) <- names(design.matrix)
results_all <- as.data.frame(cbind(design_results, results_all))
names(results_all) <- c(names(design.matrix), names(ssd_results))
save(results_all, file = "AllResults.Rdata")

# Plots
## Grouped the results by n1, n2, rho, eff.size and bf.threshold and estimate median and mean by groups.
results_all <- mutate(results_all, group = paste(n1, n2, rho, eff.size, bf.thresh, sep = ' '))
results_grouped <- results_all %>% 
    as.factor(results_all$group) %>% 
    summarise(median.BF12 = median(BF.12), median.BF21 = median(BF.21),
              mean.PMP1 = mean(PMP.1), mean.PMP2 = mean(PMP.2), mean.n2 = mean(n2), mode.n2 = mode(n2)) %>% 
    as.data.frame()
results_grouped <- cbind(design.matrix, results_grouped)

## Bayes factors -------------------------------------
hist(log10(results_grouped$median.BF12), breaks = 100)
ggplot(results_grouped, aes(median.BF12, color = eff.size, fill = eff.size)) +
    geom_histogram(alpha = 0.5, bins = 50) +
    scale_color_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
    facet_grid(rows = vars(n1), cols = vars(n2), labeller = label_both) +
    labs(title = "Bayes Factor H1 vs H2", subtitle = "H1:Dintervention>Dcontrol") +
    xlab("Bayes Factor") + ylab("Frequency") + 
    scale_x_log10(breaks = 10^(2:13), labels = trans_format("log10", math_format(10^.x)))

## sample size ---------------------------------------
ggplot(results_grouped, aes(y = mean.n2, x = n1)) + geom_line() +
    facet_grid(rows = vars(eff.size), cols = vars(rho), labeller = label_both) +
    labs(title = "Numer of clusters in function of clusters' size", 
         subtitle = "H1:Dintervention>Dcontrol") +
    xlab("Clusters' size") + ylab("Number of clusters") 
