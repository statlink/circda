\name{bic.mixcirc}
\alias{bic.mixcirc}
\title{
BIC and ICL for the model based clustering with circular distributions
}

\description{
BIC and ICL for the model based clustering with circular distributions.
}

\usage{
bic.mixcirc(u, rads = TRUE, type = "vm", G = 5, tol = 1e-4, maxiters = 500)
}

\arguments{
\item{u}{
A matrix containing directional data.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises mixtures, "cp" is circular Purkayastha mixtures,
"pn" is projected normal mixtures, "gcpc" is GCPC mixtures and "cipc" is CIPC
(or wrapped Caucy) mixtures.
}
\item{G}{
The maximum number of clusters to be tested. Default value is 5.
}
\item{tol}{
The tolerance value to terminate the EM algorithm.
}
\item{maxiters}{
The maximum number of iterations the EM algorithm will perform.
}
}

\details{
The function computes the BIC and ICL to decide on the optimal number of clusters when
using mixtures of SESPC or mixtures of ESAG distributions.
}

\value{
A plot of the ICL values and a list including:
\item{bic}{
The BIC values for all the models tested.
}
\item{icl}{
The ICL values for all the models tested.
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time,
the second element is the system time and the third element is the elapsed time.
}
}

\references{
Perdikis T., Alharbi N. and Tsagris M. (2026). Model--based clustering for spherical
and hyper--spherical data using elliptically symmetric distributions.

https://arxiv.org/abs/2605.27496
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{mixcirc.mle}, \link{circ.mle}
}
}

\examples{
u1 <- rnorm(100, 3, 0.1)
u2 <- rnorm(100, 4, 0.1)
u <- c(u1, u2)
bic.mixcirc(u, type = "vm", G = 5)
}

