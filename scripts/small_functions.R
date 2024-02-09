###### SMALL FUNCTIONS #####

# Model fitting
fit_lmer <- function(x){
    suppressMessages({
    fitted_model <- lmer(resp ~ intervention + control - 1 + (1 | id), data = x)})
    return(fitted_model)
}

# Obtain the variance covariance matrix
varcov <- function(output.lmer, number) {
    varcov <- matrix(vcov(output.lmer)[number], nrow = 1, ncol = 1)
    return(varcov)
}

# Obtain variances
get_variances <- function(output.lmer, row) {
    variances <- as.data.frame(VarCorr(output.lmer))
    value <- variances[row, 4]
    return(value)
}

# Marking singular matrices
marker_func <- function(output.lmer) {
    ifelse(isSingular(output.lmer), marker <- 1, marker <- 0)
}

# Extract results
extract_res <- function(x, number){
    results <- x[[number]]
    return(results)
}
