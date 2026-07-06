\name{rmixcirc}
\alias{rmixcirc}
\title{
Random value generation from mixtures of circular distributions
}

\description{
Random value generation from mixtures of circular distributions.
}

\usage{
rmixcirc(n, rads = TRUE, probs, mu, kappa, rho, type = "vm")
}

\arguments{
\item{n}{
The sample size.
}
\item{rads}{
The data that to be generated will be expressed in rads (TRUE) or angles (FALSE)?
}
\item{probs}{
A vector with the mixing probability of each group.
}
\item{mu}{A vector with the circular mean or location parameter (in radians) of each group.
}
\item{kappa}{
A vector with the concentration parameter of each group.
}
\item{rho}{
A vector with the \eqn{\rho} parameter of the GCPC distribution.
}
\item{type}{
The distribution to fit, "vm" is von Mises mixtures, "cp" is circular Purkayastha mixtures,
"pn" is projected normal mixtures, "gcpc" is GCPC mixtures and "cipc" is CIPC
(or wrapped Caucy) mixtures.
}
}

\details{
The function simulates random values simulation from a given mixture of circular distributions.
}

\value{
A list including:
\item{id}{
An indicator of the group of each simulated vector.
}
\item{u}{
A vector with the simulated data.
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
\code{\link{mixcirc.mle}
}
}

\examples{
x <- as.matrix( iris[, 1:3] )
x <- x / sqrt( rowSums(x^2) )
}
