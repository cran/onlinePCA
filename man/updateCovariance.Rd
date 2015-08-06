\name{updateCovariance}
\alias{updateCovariance}

\title{
Update the Sample Covariance Matrix}

\description{
This function recursively updates a covariance matrix without entirely recomputing it when new observations arrive.}

\usage{
updateCovariance(C, x, n, xbar, f, byrow = TRUE) 
}

\arguments{
  \item{C}{covariance matrix.}
  \item{x}{vector/matrix of new data.}
  \item{n}{sample size before observing \code{x}.}
  \item{xbar}{mean vector before observing \code{x}.}
  \item{f}{forgetting factor: a number beween 0 and 1.}
  \item{byrow}{Are the observation vectors in \code{x} stored in rows?}
}

\details{
The forgetting factor \code{f} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{f} must be specified. If \code{f} is specified, its value overrides the argument \code{n}. The default \code{f=1/n} corresponds to a stationnary observation process.\cr
The argument \code{byrow} should be set to TRUE (default value) if the data vectors in \code{x} are stored in rows and to FALSE if they are stored in columns. The function automatically handles the case where \code{x} is a single vector.}

\value{
The updated covariance matrix. 
}

\seealso{
\code{\link{updateMean}}
}

\examples{
n <- 1e4
n0 <- 5e3
d <- 10
mat <- matrix(runif(n*d), n, d)

## Direct computation of the covariance
C <- cov(mat)

## Recursive computation of the covariance
xbar0 <- colMeans(mat[1:n0,])
C0 <- cov(mat[1:n0,])
Crec <- updateCovariance(C0, mat[(n0+1):n,], n0, xbar0)

## Check equality
all.equal(C, Crec)
}


