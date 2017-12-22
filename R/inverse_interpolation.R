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


#' Local inverse interpolation with bandwidth selected via k-fold cross-validation.
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param ypred The matrix of "test" values of y for which a value of x will be predicted.
#' @param min_omega The minimum bandwidth attempted.
#' @param max_omega The maximum bandwidth attempted. 
#' @param k The number of folds to use in cross-validation.
#' @param offset 
#' @details The best omega is that in (min_omega, max_omega) which minimizes the
#'          median relative error among the k cross-validation tests.
#' @export
#' @return List of `omega_min`, the omega that minimizes the cross-validation error,
#'                 `ypred` the matrix of predicted x values for the optimal
#'                 omega, with dimensions (nrow(ypred), ncol(x)), and
#'                 `boundary_value` a boolean indicating if `omega_min` is on
#'                 the boundary of the range of omegas examined.
automatic_inverse_interpolation <- function (x, y, ypred, min_omega, k=5,
					     max_omega=50*min_omega, offset=0) {
    if (any(x <= 0) && any(x >= 0) && offset == 0) {
	    warning("Values of x span zero: should set nonzero offset.")
    }
    dim(x) <- c(NROW(x), NCOL(x))
    omegas <- seq(min_omega, max_omega, length.out=100)
    shuffled_idxs <- sample.int(nrow(x))
    folds <- split(shuffled_idxs, cut(seq_along(shuffled_idxs), breaks=k)) 
    median_rel_error <- sapply(omegas, function(omega) {
        rel_error <- sapply(1:length(folds), function(i) {
           test <- folds[[i]]
           train <- unname(do.call(c, folds[-i]))
    	   xpred <- inverse_interpolation(x=x[train, ], y=y[train, ], y[test, ], omega=omega)
    	   rel_err <- abs(xpred - x[test, ]) / (offset + abs(x[test, ]))
    	   sum(rel_err)
        })
	median(rel_error)
        })
    best_omega <- omegas[which.min(median_rel_error)]
    return(list(omega_min=best_omega,
		ypred=inverse_interpolation(x=x, y=y, ypred=ypred, omega=best_omega),
	        boundary_value=! (min(omegas) < best_omega && best_omega < max(omegas))
		))
}
