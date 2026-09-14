\name{circ.ridge}
\alias{circ.ridge}
\title{
Ridge circular regression
}

\description{
Ridge circular regression.
}

\usage{
circ.ridge(y, x, rads = TRUE, type = "vm", lambda = NULL, nlambda = 100,
xnew = NULL, tol = 1e-6, maxiters = 100)
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
\item{xnew}{
The new values of some independent variable(s) whose circular values you want to predict.
The can be Euclidean or categorical. If you have no new x values, leave it NULL (default).
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson algorithm will perform.
}
}

\details{
Ridge regression using the von Mises, the projected normal or the CIPC (or wrapped Cauchy) distribution.
}

\value{
A list including:
\item{runtime}{
The runtime of the ridge regression model.
}
\item{info}{
A matrix with two columns containing the lambda, and the penalized log-likelihood.
}
\item{be}{
A list with the ridge beta coefficients.
}
\item{est}{
A list with the ridge estimated compositions.
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
\code{\link{cv.circridge}, \link{circ.reg}
}
}

\examples{
y <- rcirc(100, mu = 3, kappa = 2, rads = TRUE, type = "vm")
x <- matrix( rnorm(100 * 10), ncol = 10)
circ.ridge(y, x, type = "vm")
}

