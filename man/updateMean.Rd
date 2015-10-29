\name{updateMean}
\alias{updateMean}

\title{
Update the Sample Mean Vector
}

\description{
Recursive update of the sample mean vector.
}

\usage{
updateMean(xbar, x, n, f, byrow = TRUE)
}

\arguments{
  \item{xbar}{current mean vector.}
  \item{x}{vector or matrix of new data.}
  \item{n}{sample size before observing \code{x}.}
  \item{f}{forgetting factor: a number in (0,1).}
  \item{byrow}{Are the observation vectors in \code{x} stored in rows (TRUE) or in columns (FALSE)?}
}

\details{
The forgetting factor \code{f} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{f} must be specified. If \code{f} is specified, its value overrides the argument \code{n}. For a given argument \code{n}, the default value of\code{f} is\eqn{k/(n+k)}, with \eqn{k} the number of new vector observations. This corresponds to a stationnary observation process.\cr
% The argument \code{byrow} should be set to TRUE (default value) if the data vectors in \code{x} are stored in rows and to FALSE if they are stored in columns. The function automatically handles the case where \code{x} is a single vector.
}

\value{
The updated mean vector. 
}

\seealso{
\code{\link{updateCovariance}}
}

\examples{
n <- 1e4
n0 <- 5e3
d <- 10
x <- matrix(runif(n*d), n, d)

## Direct computation
xbar1 <- colMeans(x)

## Recursive computation
xbar2 <- colMeans(x[1:n0,])
xbar2 <- updateMean(xbar2, x[(n0+1):n,], n0)

## Check equality
all.equal(xbar1, xbar2)
}

