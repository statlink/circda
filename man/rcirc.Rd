\name{rcirc}
\alias{rcirc}
\title{
Random value generation from circular distributions
}

\description{
Random value generation from circular distributions.
}

\usage{
rcirc(n, rads = TRUE, mu, kappa, rho, type = "vm")
}

\arguments{
\item{n}{
The sample size.
}
\item{rads}{
The data that to be generated will be expressed in rads (TRUE) or angles (FALSE)?
}
\item{mu}{A vector with the circular mean or location parameter (in radians).
}
\item{kappa}{
A vector with the concentration parameter.
}
\item{rho}{
A vector with the \eqn{\rho} parameter of the GCPC distribution.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, "gcpc" is GCPC distribution and "cipc" is CIPC
(or wrapped Caucy) distribution.
}
}

\details{
The function simulates random values simulation from a given circular distribution.
}

\value{
A vector with the simulated data.
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
\code{\link{circ.mle}, \link{rmixcirc}
}
}

\examples{
x <- as.matrix( iris[, 1:3] )
x <- x / sqrt( rowSums(x^2) )
}
