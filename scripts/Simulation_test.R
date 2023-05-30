############################### TEST #########################################

# Design matrix
n1 <- c(5, 10, 20, 40)
n2 <- c(30, 60 ,90)
rho <- c(0.25, 0.05, 0.1)
eff.size <- c(0.2, 0.5, 0.8)
design.matrix <- expand.grid(n1, n2, rho, eff.size)

nrow.design <- nrow(design.matrix)
iterations <- 10

for (Row in seq(nrow.design)) {
    for (i in seq(iterations)) {
        
    }
}




# - Add plots.This could be a function.