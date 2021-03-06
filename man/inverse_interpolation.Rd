% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inverse_interpolation.R
\name{inverse_interpolation}
\alias{inverse_interpolation}
\title{Local inverse interpolation}
\usage{
inverse_interpolation(x, y, ypred, omega)
}
\arguments{
\item{x}{The matrix of "training" values of x, with one row per training set.}

\item{y}{The matrix of "training" values of y corresponding to x.}

\item{ypred}{The matrix of "test" values of y for which a value of x will be predicted.}

\item{omega}{The bandwidth on the local weights (can be a vector of length equal to the number of columns of y).}
}
\value{
The matrix of predicted x values, with dimensions (nrow(ypred), ncol(x)).
}
\description{
Local inverse interpolation
}
