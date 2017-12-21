#' Local inverse interpolation
#'
#' @param x The matrix of "training" values of x, with one item per training set.
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


#' Local inverse interpolation with bandwidth selected via k-fold cross-validation.
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param ypred The matrix of "test" values of y for which a value of x will be predicted.
#' @param seed The random seed used in the k-fold cross-validation.
#' @param min_omega The minimum bandwidth attempted. 
#' @param max_omega The initial maximum bandwidth attempted. 
#' @param k The number of folds to use in cross-validation.
#' @export
#' @return List of `omega_min`, the omega that minimizes the cross-validation error,
#'                 `ypred` the matrix of predicted x values for the optimal
#'                 omega, with dimensions (nrow(ypred), ncol(x)), and
#'                 `boundary_value` a boolean indicating if `omega_min` is on
#'                 the boundary of the range of omegas examined.
automatic_inverse_interpolation <- function (x, y, ypred, seed, min_omega, max_omega=NA, k=5) {
    max_omega <- 50 * min_omega
    omegas <- seq(min_omega, max_omega, length.out=100)
    idxs <- 1:length(x)
    set.seed(seed)
    shuffled_idxs <- sample(idxs, length(idxs), replace=FALSE)
    folds <- split(shuffled_idxs, cut(shuffled_idxs, breaks=k))
    best_omegas <- sapply(1:length(folds), function(i) {
        test <- folds[[i]]
        train <- unname(do.call(c, folds[-i]))
    	errors <- sapply(omegas, function(omega) {
    	    xpred <- inverse_interpolation(x=x[train], y=y[train, ], y[test, ], omega=omega)
    	    rel_err <- abs(xpred - x[test]) / x[test]
    	    sum(rel_err)
            })
    	best_omega <- omegas[which(min(errors) == errors)]
    	best_omega
        })
    best_omega <- median(do.call(c, as.list(best_omegas)))
    return(list(omega_min=best_omega,
		ypred=inverse_interpolation(x=x, y=y, ypred=ypred, omega=best_omega),
	        boundary_value=! (min(omegas) < best_omega && best_omega < max(omegas))
		))
}
