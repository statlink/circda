\name{von Mises-Fisher kernel density estimation for (hyper-)spherical data}
\alias{vmf.kernel}
\title{
Kernel density estimation for (hyper-)spherical data using a von Mises-Fisher kernel
}

\description{
A von Mises-Fisher kernel is used for the non parametric density estimation.
}

\usage{
vmf.kernel(y, h = "rot", ncores = 1)
}

\arguments{
\item{y}{
A matrix with unit vectors, i.e. the data being expressed in Euclidean coordinates.
}
\item{h}{
The bandwidth. Either "rot" for the rule suggested by Garcia-Portugues (2013) or a number.
}
\item{ncores}{
The number of cores (threads) to use. The default is 1. If the package was built without
OpenMP support, this argument is ignored and a single core is used.
}
}

\details{
A von Mises-Fisher kernel is used for the non parametric density estimation.
The density estimate at each observation requires all pairwise inner products, which are
computed in C++. When ncores > 1 the computation is split across threads using OpenMP.
}

\value{
A list including:
\item{h}{
The bandwidth used.
}
\item{f}{
A vector with the kernel density estimate calculated for each of the data points.
}
}

\references{
Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for
kernel density estimation with directional data.
Electronic Journal of Statistics, 7, 1655--1685.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vm.kernel}
}
}

\examples{
y <- Rfast::rvmf(1000, rnorm(5), 15)
f <- vmf.kernel(y, h = "rot")
f2 <- vmf.kernel(y, h = "rot", ncores = 2)
}
