######################### Practtical example in paper ##########################

var_e <- 45
var_u <- 3.5
delta <- 1.39 # Unstandardised effect size
n1 <- 30

# Number of clusters based on ICC and standardised effect size
ICC <- (var_u) / (var_e + var_u)
delta_st <- 1.39 / sqrt(var_e + var_u) # standardised effect size
n2 <- 4*((1 + (n1 - 1) * ICC) / n1) * ((qnorm(0.975) + qnorm(0.8)) / delta_st)^2

# Equation from: 
#   Moerbeek, M., & Teerenstra, S. (2016). Power analysis of trials with multilevel
#       data. CRC Press. http://www.crcnetbase.com/isbn/9781498729901
# Example from:
#   Moerbeek M. (2006). Power and money in cluster randomized trials: when is it
#       worth measuring a covariate?. Statistics in medicine, 25(15), 2607–2617.
#       https://doi.org/10.1002/sim.2297