\name{vm.kernel.da}
\alias{vm.kernel.da}
\title{
Kernel discriminant analysis of circular data with a von Mises kernel
}

\description{
Kernel discriminant analysis of circular data with a von Mises kernel.
}

\usage{
vm.kernel.da(ynew, y, ina, h = "tay", rads = TRUE, ncores = 1)
}

\arguments{
\item{ynew}{
A numeric vector containing the new circular data.
}
\item{y}{
A numeric vector containing the circular data.
}
\item{ina}{
A numeric vector (or factor) with the group of each observation in y, taking the values 1, 2, ...
}
\item{h}{
The bandwidth. Either "tay" for the rule suggested by Taylor (2008), "rot" for the rule
suggested by Garcia-Portugues (2013), a single number used for all groups, or a vector with
one bandwidth for each group.
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
(cross-validation) or estimate it as Taylor (2008) or Garcia-Portugues (2013), separately for each group.

Each new observation is allocated to the group with the highest score, i.e. the logarithm of
the prior probability of the group times the kernel density estimate of that group at the new
observation. The scores are computed in C++, one new observation at a time, so the matrix of
all inner products between the new and the training observations is never formed and the memory
requirements stay low even when there are many new observations. When ncores > 1 the new observations are split across threads using OpenMP.
}

\value{
A list including:
\item{h}{
The estimated bandwidth of each group.
}
\item{class}{
The estimated class of the new data.
}
\item{score}{
A matrix with the score of each new observation (rows) for each group (columns).
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
\code{\link{vm.kernel}, \link{vmf.kernel.da}, \link{vm.nb}
}
}

\examples{
y1 <- Rfast::rvonmises(100, 1, 5, rads = TRUE)
y2 <- Rfast::rvonmises(100, 3, 5, rads = TRUE)
y <- c(y1, y2)
ina <- rep(1:2, each = 100)
ynew <- Rfast::rvonmises(50, 1, 5, rads = TRUE)
mod <- vm.kernel.da(ynew, y, ina, h = "tay")
table(mod$class)
mod2 <- vm.kernel.da(ynew, y, ina)
}
