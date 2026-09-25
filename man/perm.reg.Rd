\name{perm.reg}
\alias{perm.reg}
\title{
Permutation based log-likelihood ratio test for one or more predictor variables
}

\description{
Permutation based log-likelihood ratio test for one or more predictor variables
}

\usage{
perm.reg(y, x, rads = TRUE, type = "vm", tol = 1e-6, maxiters = 500, R = 999)
}

\arguments{
\item{y}{
A vector with the circular data expressed in radians, or angles.
}
\item{x}{
A numerical matrix with many variables. A circular regression will be fit to each of these variables.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, and "cipc" is CIPC (or wrapped Cauchy) distribution.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson algorithm will perform.
}
\item{R}{
The number of permutations to perform.
}
}

\details{
The function performs the permutation-based log-likelihood ratio test to test the significance
of one or more predictor variables against the model with no predictors.
}

\value{
A vectror with the test statistic and the p-value.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. Journal of the American Statistical Association, 93(443): 1068--1077.

Tsagris M. and Alzeley O. (2025). Circular and spherical projected Cauchy distributions:
A Novel Framework for Circular and Directional Data Modelling.
Australian & New Zealand Journal of Statistics, 67(1): 77--103.
https://arxiv.org/pdf/2302.02468.pdf

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Max-
imum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.reg}
}
}

\examples{
y <- rcirc(100, mu = 3, kappa = 2, rads = TRUE, type = "vm")
x <- rnorm(100)
perm.reg(y, x, rads = TRUE, type = "vm", R = 499)
}
