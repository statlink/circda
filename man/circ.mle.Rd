\name{circ.mle}
\alias{circ.mle}
\title{
MLE of a circular distribution
}

\description{
MLE of a circular distribution.
}

\usage{
circ.mle(u, rads = TRUE, type ="vm", tol = 1e-6, maxiters = 100)
}

\arguments{
\item{u}{
A vector with the circular data expressed in radians, or angles.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, "gcpc" is GCPC distribution and "cipc" is CIPC
(or wrapped Caucy) distribution.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson algorithm will perform.
}
}

\details{
The functions performs model-based clustering using mixtures of circular distributions.
The initial step of the algorithm is based on the k-means algorithm.
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
\code{\link{bic.mixcirc}
}
}

\examples{
x <- as.matrix( iris[, 1:3] )
x <- x / sqrt( rowSums(x^2) )
}
