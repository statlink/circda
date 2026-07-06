\name{dmixcirc}
\alias{dmixcirc}
\title{
Density values for mixtures of circular distributions
}

\description{
Density values for mixtures of circular distributions.
}

\usage{
dmixcirc(u, rads = TRUE, probs, mu, kappa, rho, type = "vm", logden = FALSE)
}

\arguments{
\item{u}{
A vector with the values at which the density is to be computed.
}
\item{rads}{
The data that to be generated will be expressed in rads (TRUE) or angles (FALSE)?
}
\item{probs}{
A vector with the mixing probability of each group.
}
\item{mu}{A vector with the circular mean (in radians) of each group.
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
\item{logden}{
If you want the log of the density set this equal to TRUE.
}
}

\details{
The function computes the density of a mixture of circular distributions at the given points.
}

\value{
A vector with the density values of a mixture of circular distributions computed at the given points.
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
