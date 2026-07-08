\name{circ.reg}
\alias{circ.reg}
\title{
Regression model using a circular distribution
}

\description{
Regression model using a circular distribution.
}

\usage{
circ.reg(y, x, rads = TRUE, type = "vm", xnew = NULL, tol = 1e-6, maxiters = 100)
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
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, "gcpc" is GCPC distribution and "cipc" is CIPC
(or wrapped Caucy) distribution.
}
\item{xnew}{
The new values of some independent variable(s) whose circular values you want to predict.
The can be Euclidean or categorical. If you have no new x values, leave it NULL (default).
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations allowed in the Newton-Raphson algorithm.
}
}

\details{
The functions performs regression using circular distributions.
}

\value{
A list including:
\item{param}{
A matrix with the mixing probability of each group and the estimated parameters of
the chosen distribution.
}
\item{loglik}{
The value of the maximised log-likelihood of the chosen distribution.
}
\item{iter}{
The number of iteration required by the EM algorithm.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Perdikis T., Alharbi N. and Tsagris M. (2026). Model--based clustering for spherical
and hyper--spherical data using elliptically symmetric distributions.

https://arxiv.org/abs/2605.27496
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.mle}
}
}

\examples{
x <- as.matrix( iris[, 1:3] )
x <- x / sqrt( rowSums(x^2) )
}
