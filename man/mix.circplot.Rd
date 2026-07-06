\name{mix.circplot}
\alias{mix.circplot}
\title{
Circular histogram
}

\description{
Histogram of circular data.
}

\usage{
mix.circplot(u, ina, mu, rads = TRUE, lwd = 2)
}

\arguments{
\item{u}{
A vector with the circular data expressed in radians, or angles.
}
\item{ina}{
A numeric vector indicating the groups of the data.
}
\item{mu}{
The estimated circular mean or a circular location estimate. If this
is not available the circular mean is caclulated internally.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{lwd}{
The thickness of the above arrow of the circular means.
}
}

%\details{
%}

\value{
A circular histogram ith different colours for each group and each circular mean.
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{circplot}
}
}

\examples{
u <- c( rnorm(50, rnorm(3, 0.2) ), rnorm(50, rnorm(1, 0.2) ) )
ina <- rep(1:2, each = 50)
mix.circplot(u, ina = ina, mu = c(3, 1) )
}

