\name{colcirc.mle}
\alias{colcirc.mle}
\title{
Column-wise MLE of some circular distributions
}

\description{
Column-wise MLE of some circular distributions.
}

\usage{
colcirc.mle(u, rads = TRUE, type = "vm", tol = 1e-07, maxiters = 100, parallel = FALSE)
}

\arguments{
\item{u}{
A numerical matrix with the circular data.
}
\item{rads}{
If the data are expressed in angles set this to FALSE.
}
\item{type}{
The distribution to fit, "vm" is von Mises distribution, "cp" is the circular Purkayastha distribution,
"pn" is projected normal distribution, "gcpc" is GCPC distribution and "cipc" is CIPC
(or wrapped Caucy) distribution.
}
\item{tol}{
The tolerance level to stop the iterative process of finding the MLEs.
}
\item{maxiters}{
The maximum number of iterations to implement. This is for the "spml" only.
}
\item{parallel}{
Should the computations take place in parallel? This is for the "spml" only.
}
}

\details{
The function performs MLE of circular distributions for each column of a matrix.
}

\value{
A matrix with two, columns. The first one contains the parameters of the distribution
and the second columns contains the log-likelihood values.
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of \eqn{I_s(x)}. Computational Statistics, 27(1): 177--190.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. Journal of the American Statistical Association, 93(443): 1068--1077.

Tsagris M. and Alzeley O. (2025). Circular and spherical projected Cauchy distributions:
A Novel Framework for Circular and Directional Data Modelling.
Australian & New Zealand Journal of Statistics, 67(1): 77--103.
https://arxiv.org/pdf/2302.02468.pdf

Alzeley O. and Tsagris M. (2026). On the generalized circular projected Cauchy distribution.
Mathematics, 14(11): 1934. https://www.mdpi.com/2227-7390/14/11/1934

Purkayastha S. (1991). A Rotationally Symmetric Directional Distribution: Obtained through Max-
imum Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70--83.

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. Communications
in Statistics-Theory and Methods, 19(6): 1973--1986
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{circ.mle}
}
}

\examples{
u <- matrix( circda::rcirc(200 * 10, mu = 3, kappa = 5), ncol = 10 )
res <- colcirc.mle(u)
}

