#' Local inverse interpolation
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param ypred The matrix of "test" values of y for which a value of x will be predicted.
#' @param omega The bandwidth on the local weights (can be a vector of length equal to the number of columns of y).
#' @export
#' @return The matrix of predicted x values, with dimensions (nrow(ypred), ncol(x)).
inverse_interpolation <- function (x, y, ypred, omega) {
    stopifnot(NCOL(y) == NCOL(ypred))
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


#' Local inverse interpolation with bandwidth that minimizes relative error
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param min_omega The minimum bandwidth attempted.
#' @param max_omega The maximum bandwidth attempted. 
#' @param k The number of folds to use in cross-validation.
#' @param offset 
#' @details The best omega is that in (min_omega, max_omega) which minimizes the
#'          median relative error among the k cross-validation tests.
#' @export
#' @return List of `omega_min`, the omega that minimizes the cross-validation error,
#'                 `error_min`, the minimum relative error,
#'                 `predict` a function of `ypred` that will call
#'                 inverse_interpolation to return the matrix of predicted x
#'                 values for the optimal omega, with dimensions (nrow(ypred),
#'                 ncol(x)), and `boundary_value` a boolean indicating if
#'                 `omega_min` is on the boundary of the range of omegas
#'                 examined.
inverse_interpolation_minRE <- function (x, y, min_omega, k=5,
					     max_omega=50*min_omega, offset=0) {
    if (any(x <= 0) && any(x >= 0) && offset == 0) {
	    warning("Values of x span zero: should set nonzero offset.")
    }
    automatic_inverse_interpolation_general(x, y, min_omega,
					    loss=function(yh, y) abs(yh - y) / (offset + abs(y))
					    ,
					    k=k, max_omega=max_omega)
}

#' Local inverse interpolation with bandwidth selected via k-fold cross-validation.
#'
#' @param x The matrix of "training" values of x, with one row per training set.
#' @param y The matrix of "training" values of y corresponding to x.
#' @param min_omega The minimum bandwidth attempted.
#' @param max_omega The maximum bandwidth attempted. 
#' @param k The number of folds to use in cross-validation.
#' @param loss A function(yhat, y) the gives the error.
#' @details The best omega is that in (min_omega, max_omega) which minimizes the
#'          median sum of the loss among the k cross-validation tests.
#' @export
#' @return List of `omega_min`, the omega that minimizes the cross-validation error,
#'                 `error_min`, the minimum relative error,
#'                 `predict` a function of `ypred` that will call
#'                 inverse_interpolation to return the matrix of predicted x
#'                 values for the optimal omega, with dimensions (nrow(ypred),
#'                 ncol(x)), and `boundary_value` a boolean indicating if
#'                 `omega_min` is on the boundary of the range of omegas
#'                 examined.
automatic_inverse_interpolation_general <- function (x, y, min_omega, loss, k=5, max_omega=50*min_omega) {
    dim(x) <- c(NROW(x), NCOL(x))
    omegas <- seq(min_omega, max_omega, length.out=100)
    shuffled_idxs <- sample.int(nrow(x))
    folds <- split(shuffled_idxs, cut(seq_along(shuffled_idxs), breaks=k)) 
    median_loss <- sapply(omegas, function(omega) {
        total_loss <- sapply(1:length(folds), function(i) {
           test <- folds[[i]]
           train <- unname(do.call(c, folds[-i]))
    	   xpred <- inverse_interpolation(x=x[train, ], y=y[train, ], y[test, ], omega=omega)
    	   sum(loss(xpred, x[test, ]))
        })
	median(total_loss)
        })
    best_omega <- omegas[which.min(median_loss)]
    return(list(omega_min=best_omega,
		error_min=min(median_loss),
		predict=function(ypred) inverse_interpolation(x=x, y=y, ypred=ypred, omega=best_omega),
	        boundary_value=! (min(omegas) < best_omega && best_omega < max(omegas))
		))
}
