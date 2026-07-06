\name{circ.da}
\alias{circ.da}
\title{
Prediction of a new observation using discriminant analysis based on circular distributions
}

\description{
Prediction of a new observation using discriminant analysis based on circular distributions.
}

\usage{
circ.da(unew, u, ina, rads = TRUE, type = c("vm", "cp", "pn", "gcpc", "cipc") )
}

\arguments{
\item{unew}{
The new observation(s) (expressed in radians, or angles) whose group is to be predicted.
}
\item{u}{
A numerical vector with the data (expressed in radians, or angles).
}
\item{ina}{
A vector indicating the groups of the data y.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, "gcpc" is GCPC distribution and "cipc" is CIPC
(or wrapped Caucy) distribution.
}
}

\details{
Prediction of the class of a new circular vector assuming some distributions.
}

\value{
A matrix with 5 columns, one for each distribution, where each row contains the predicted class
of each observation If you chose less distributions to test, the columns of the non-chosen
distributions will contain NAs.
}

\references{
Tsagris M. and Alenazi A. (2019). Comparison of discriminant analysis methods on the sphere.
Communications in Statistics: Case Studies, Data Analysis and Applications, 5(4): 467--491.

Tsagris M., Papastamoulis P. and Kato S. (2025). Directional data analysis using the spherical Cauchy and the Poisson kernel-based distribution.
Statistics and Computing, 35:51.

Morris J. E. and Laycock P. J. (1974). Discriminant analysis of directional data. Biometrika, 61(2): 335--341.

Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{\link{circda.cv}
}
}

\examples{
m1 <- rnorm(3, 0.1)
m2 <- rnorm(5, 0.1)
}
