#####################################################
# Formatting formula for p value
format.p <- function(p, precision = .001) {
    digits <- -log(precision, base = 10)
    p <- formatC(p, format = 'f', digits = digits)
    if (p < .001) {
        p = paste0('< ', precision)}
    if (p >= .001) {
        p = paste0('= ', p)    }
    sub("", "", p)
}
######################################################