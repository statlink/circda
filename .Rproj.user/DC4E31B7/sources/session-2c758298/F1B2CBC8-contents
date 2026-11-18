\name{von Mises kernel density estimation}
\alias{vm.kernel}
\title{
Kernel density estimation of circular data with a von Mises kernel
}

\description{
Kernel density estimation of circular data with a von Mises kernel.
}

\usage{
vm.kernel(y, h = "tay", rads = TRUE, ncores = 1)
}

\arguments{
\item{y}{
A numeric vector containing the circular data.
}
\item{h}{
The bandwidth. Either "tay" for the rule suggested by Taylor (2008), "rot" for the rule
suggested by Garcia-Portugues (2013) or a number.
}
\item{rads}{
If the data are in radians, this should be TRUE and FALSE otherwise.
}
\item{ncores}{
The number of cores (threads) to use. The default is 1. If the package was built without
OpenMP support, this argument is ignored and a single core is used.
}
}

\details{
The user has the option to use a bandwidth he/she has found in some way
(cross-validation) or estimate it as Taylor (2008) or Garcia-Portugues (2013).

The density estimate at each observation requires all pairwise inner products, which are
computed in C++. When ncores > 1 the computation is split across threads using OpenMP.
}

\value{
A list including:
\item{h}{
The bandwidth. If the user chose one of "tay" or "rot" the estimated bandwidth will be returned.
}
\item{f}{
The kernel density estimate at the observed points.
}
}

\references{
Taylor, C. C. (2008). Automatic bandwidth selection for circular density estimation.
Computational Statistics & Data Analysis, 52(7): 3493--3500.

Garcia Portugues, E. (2013). Exact risk improvement of bandwidth selectors for kernel
density estimation with directional data. Electronic Journal of Statistics, 7, 1655--1685.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{vmf.kernel}
}
}

\examples{
y <- Rfast::rvonmises(1000, 2.4, 10, rads = TRUE)
hist(y, freq = FALSE)
f1 <- vm.kernel(y, h = "rot", rads = TRUE)$f
f2 <- vm.kernel(y, h = "tay", rads = TRUE)$f
points(y, f1, col = 1)
points(y, f2, col = 2)
}
