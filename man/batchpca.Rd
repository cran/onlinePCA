\name{batchpca}
\alias{batchpca}

\title{Batch PCA}

\description{
This function performs the PCA of a data matrix or covariance matrix, returning the specified number of principal components (eigenvectors) and eigenvalues.}

\usage{
batchpca(x, q, center, type = c("data","covariance"), byrow = FALSE)
}

\arguments{
  \item{x}{data or covariance matrix}
  \item{q}{number of requested PCs}
  \item{center}{optional centering vector for \code{x}}
  \item{type}{type of the matrix \code{x}}
  \item{byrow}{Are observation vectors stored in rows (TRUE) or in columns (FALSE)?}
}

\details{
The PCA is efficiently computed using the functions \code{\link[RSpectra]{svds}} or \code{\link[RSpectra]{eigs_sym}} of package \code{RSpectra}, depending on the argument \code{type}. An Implicitly Restarted Arnoldi Method (IRAM) is used in the former case and an Implicitly Restarted Lanczos Method (IRLM) in the latter.\cr
The arguments \code{center} and \code{byrow} are only in effect if \code{type} is \code{"data"}. In this case a scaling factor \eqn{1/\sqrt{n}} (not \eqn{1/\sqrt{n-1}}) 
is applied to \code{x} before computing its singular values and vectors, where \eqn{n} is the number of observation vectors stored in \code{x}.    
}

\value{
A list with components
\item{values}{the first \code{q} squared singular values of \code{x} if \code{type="data"}; the first \code{Q} eigenvalues if \code{type="covariance"}.}
\item{vectors}{the first \code{q} PC of \code{x}.}
}

\references{
\url{https://www.arpack.org}
}

\examples{

\dontrun{
## Simulate data
n <- 1e4
d <- 500
q <- 10
x <- matrix(runif(n*d), n, d)
x <- x \%*\% diag(sqrt(12*(1:d)))
# The eigenvalues of cov(x) are approximately 1, 2, ..., d
# and the corresponding eigenvectors are approximately  
# the canonical basis of R^p

## PCA computation (from fastest to slowest)
system.time(pca1 <- batchpca(scale(x,scale=FALSE), q, byrow=TRUE))
system.time(pca2 <- batchpca(cov(x), q, type="covariance"))
system.time(pca3 <- eigen(cov(x),TRUE))
system.time(pca4 <- svd(scale(x/sqrt(n-1),scale=FALSE), 0, q))
system.time(pca5 <- prcomp(x))
}
}
