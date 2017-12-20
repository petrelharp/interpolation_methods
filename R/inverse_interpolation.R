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
    ypred <- sweep(ypred, 2, omega, "/")
    dvec <- fields::rdist(y, ypred)
    dvec <- sweep(dvec, 2, apply(dvec,2,min), "-")
    w <- exp(-dvec^2)
    wsums <- colSums(w)
    if (any(wsums == 0)) {
        warning("Bandwidth is too small -- some results will be NaN.")
    }
    w <- sweep(w, 2, wsums, "/")
    xpred <- crossprod(w, x)
    return(xpred)
}

