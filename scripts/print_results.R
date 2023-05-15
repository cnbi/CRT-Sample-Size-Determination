############################ OUTPUT RESULTS ##################################
# Set class name and the type
setClass("SSD")

print.SSD <- function(object) {
    setMethod("show", "SSD", function(object){
        title <- "Final sample size"
        cat(paste("\n", title, "\n", sep = ""))
        row <- paste(rep("=", nchar(title)), collapse = "")
        cat(row, "\n")
        cat("Using cluster size = ", object$n1, " and number of clusters = ", object$n2, "\n")
        cat("P (BF.12 > BF.threshold | H.1 = ", object$prop.BF12)
    })
}
# How the output would be
# setMethod("show", "SSD", function(object){
#     title <- "Final sample size"
#     cat(paste("\n", title, "\n", sep = ""))
#     row <- paste(rep("=", nchar(title)), collapse = "")
#     cat(row, "\n")
#     cat("Using cluster size = ", object$n1, " and number of clusters = ", object$n2, "\n")
#     cat("P (BF.12 > BF.threshold | H.1 = ", object$prop.BF12)
# })
