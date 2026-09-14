\name{cv.circridge}
\alias{cv.circridge}
\title{
Cross-validation for the ridge circular regression
}

\description{
Cross-validation for the ridge circular regression.
}

\usage{
cv.circridge(y, x, rads = TRUE, type = "vm", lambda = NULL, nlambda = 100, tol = 1e-06,
maxiters = 100, folds = NULL, nfolds = 10, seed = NULL)
}

\arguments{
\item{y}{
A vector with the circular data expressed in radians, or angles.
}
\item{x}{
The independent variable(s). Can be Euclidean or categorical (factor variables).
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "pn" is projected normal distribution,
"cipc" is CIPC (or wrapped Cauchy) distribution.
}
\item{lambda}{
A vector with a sequence of the ridge \eqn{\lambda} values. If you do not know it leave it NULL
and the function will compute it internally.
}
\item{nlambda}{
The number of \eqn{\lambda} values to compute.
}
\item{tol}{
The tolerance value to terminate the Fisher scoring algorithm.
}
\item{maxiters}{
The maximum number of iterations the Fisher scoring algorithm will perform.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{nfolds}{
The number of folds to produce.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
K-fold cross-validation for the ridge circular regression models.
}

\value{
A list including:
\item{runtime}{
The runtime of the ridge regression model.
}
\item{info}{
A matrix with two columns containing the lambda and the fit, the \eqn{\sum_i y_i^T\hat{y}_i/n_k}.
}
\item{lambda_opt}{
The optimal value of the ridge parameter, \eqn{\lambda}.
}
\item{be}{
The coefficients of the ridge circular regression model corresponding to the optimal \eqn{\lambda} value.
}
}

\references{
Tsagris M., Papastamoulis P. and Kato S. (2025). Directional data analysis using the spherical Cauchy and the
Poisson kernel-based distribution.
Statistics and Computing, 35:51.

Presnell B., Morrison S. P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068--1077.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{circ.ridge}
}
}

\examples{
y <- rcirc(100, mu = 3, kappa = 2, rads = TRUE, type = "vm")
x <- matrix( rnorm(100 * 10), ncol = 10)
circ.ridge(y, x, type = "vm")
}

