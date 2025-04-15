############ Script for figure in manuscript ###############

# Libraries
library(ggplot2)
library(cowplot) #Side-by-side plot

# Prior parameters
prior_mean <- 0             # Null hypothesis
prior_sd <- sqrt(1)

# Sample parameters
data_mean <- 1.5
data_sd <- 0.8

null_value <- 0
eta <- seq(-3, 3, length.out = 1000)
prior <- function(x) dnorm(x, mean = prior_mean, sd = prior_sd) # function for prior
likelihood <- function(x) dnorm(x, mean = data_mean, sd = data_sd) #likelihoood function
posterior <- function(x) prior(x) * likelihood(x)   # Baye's theorem
posterior_normalized <- posterior(eta) / sum(posterior(eta) * diff(eta)[1])

# Data with priors, likelihood and posteriors
df <- data.frame(
    treatment_eff = eta,
    Prior = prior(eta),
    Likelihood = likelihood(eta),
    Posterior = posterior_normalized
)

# Density at 0
density_0 <- df[round(df[ , 'treatment_eff'], 2) == 0.00, "Prior"]
density_0 <- (density_0[1] + density_0[2])/2
max_density_y_label <- substr(as.character(density_0), 1, 5)

# Prior plot
prior_plot <- ggplot(df, aes(x = treatment_eff)) +
    geom_line(aes(y = Prior)) +
    # Vertical line
    geom_vline(xintercept = null_value, linetype = "dotted", color = "black", size = 1) +
    xlab("Treatment effect") + ylab("Prior density") +
    theme_classic() +
    # Shaded area
    geom_ribbon(data = subset(df, treatment_eff > 0), 
                aes(ymin = 0, ymax = Prior), fill = "grey", alpha = 0.5) + 
    # Mark maximum density
    scale_y_continuous(breaks = c(0.0, 0.1, 0.2, 0.3, as.numeric(max_density_y_label), seq(0.5, 0.7, by = 0.1)),
                       limits = c(0, 0.7)) +
    scale_x_continuous(breaks = c(-3:3)) +
    # Horizontal line
    geom_hline(yintercept = max_density_y, linetype = "dotted", color = "black", size = 1)

# Density at 0
density_0 <- df[round(df[ , 'treatment_eff'], 2) == 0.00, "Posterior"]
density_0 <- (density_0[1] + density_0[2])/2
max_density_y_label <- substr(as.character(density_0), 1, 5)

# Posterior plot
posterior_plot <- ggplot(df, aes(x = treatment_eff)) +
    geom_line(aes(y = Posterior)) +
    # Vertical line
    geom_vline(xintercept = null_value, linetype = "dotted", color = "black", linewidth = 1) +
    xlab("Treatment effect") + ylab("Posterior density") +
    theme_classic() + 
    # Shaded area
    geom_ribbon(data = subset(df, treatment_eff > 0), 
                aes(ymin = 0, ymax = Posterior), fill = "grey", alpha = 0.5) + 
    # Mark density at 0 in y-axis
    scale_y_continuous(breaks = c(0.0, 0.1, as.numeric(max_density_y_label), seq(0.3, 0.8, by=0.1)),
                       limits = c(0, 0.7)) +
    scale_x_continuous(breaks = c(-3:3)) +
    # Horizontal line
    geom_hline(yintercept = density_0, linetype = "dotted", color = "black", size = 1)

# Plots side by side
plot_grid(prior_plot, posterior_plot)
ggsave("plots_fit_comp.png", width = 6, height = 4, dpi = 600)

# Bayes factor
1 - pnorm(0, mean = data_mean, sd = data_sd) #Fit H1
(1 - pnorm(0, mean = data_mean, sd = data_sd))/0.5 #BF1u
dnorm(0) # Complexity H0
density_0 # Fit H0
density_0/dnorm(0) #BF0u

((1 - pnorm(0, mean = data_mean, sd = data_sd))/0.5)/(density_0/dnorm(0)) #BF10