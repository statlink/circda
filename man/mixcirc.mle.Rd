\name{mixcirc.mle}
\alias{mixcirc.mle}
\title{
Fitting mixtures of circular distributions
}

\description{
Fitting mixtures of circular distributions.
}

\usage{
mixcirc.mle(u, type ="vm", rads = TRUE, g = 2, tol = 1e-4, maxiters = 500)
}

\arguments{
\item{u}{
A vector with the circular data expressed in radians, or angles.
}
\item{type}{
The distribution to fit, "vm" is von Mises mixtures, "cp" is circular Purkayastha mixtures,
"pn" is projected normal mixtures, "gcpc" is GCPC mixtures and "cipc" is CIPC
(or wrapped Caucy) mixtures.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{g}{
The number of groups to fit. It must be greater than or equal to 2.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxiters}{
The maximum number of iterations the EM algorithm will perform.
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
\item{probs}{
The estimated probabilities of each observation to belong to each cluster.
}
\item{pred}{
The predicted group of each observation.
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
u1 <- runif(100)
u2 <- runif(100, 0.5, 1.5)
u <- c(u1, u2)
mixcirc.mle(u, type = "vm")

}
