#' Local inverse interpolation
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param ypred The matrix of "test" values of y for which a value of x will be predicted.
#' @param omega The bandwidth on the local weights (can be a vector of length equal to the number of columns of y).
#' @export
#' @return The matrix of predicted x values, with dimensions (nrow(ypred), ncol(x)).
inverse_interpolation <- function (x, y, ypred, omega) {
    y <- sweep(y, 2, omega, "/")
    ypred <- sweep(y, 2, omega, "/")
    dvec <- fields::rdist(y, ypred)
    w <- exp(-dvec^2)
    w <- sweep(w, 2, colSums(w), "/")
    xpred <- crossprod(w, x)
    return(xpred)
}

