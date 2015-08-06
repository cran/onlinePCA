\name{batchpca}
\alias{batchpca}

\title{
Batch PCA of a covariance matrix}

\description{
This function performs the PCA of a covariance matrix, returning an arbitrary number of principal components (eigenvectors) and corresponding eigenvalues.}

\usage{
batchpca(C, q)
}

\arguments{
  \item{C}{covariance matrix}
  \item{q}{number of PCs to compute}
}

\details{
If \code{q} is small (\eqn{k <= ncol(C)/10}), the eigendecomposition is efficiently carried out by the function \code{\link[rARPACK]{eigs_sym}} of package \code{rARPACK} using an Implicitly Restarted Lanczos Method (IRLM). Otherwise the eigendecomposition is performed with the R base function \code{\link[base]{eigen}}.}

\value{
A list with components
\item{values}{the first \code{q} eigenvalues of \code{C}}
\item{vectors}{the first \code{q} PC/eigenvectors of \code{C}}
}

\references{
\url{http://www.caam.rice.edu/software/ARPACK/}
}

\examples{
n <- 1e4
p <- 500
q <- 10

X <- matrix(runif(n*p), n, p)
X <- X * rep(sqrt(12 * (1:p)), each = n)
C <- cov(X)
# The eigenvalues of cov(X) should be close to 1, 2, ..., p
# and the corresponding eigenvectors should be close to 
# the canonical basis of R^p

## q first eigenpairs of C using base function eigen
system.time(replicate(10, {
  eigenC1 <<- eigen(C, TRUE)
  eigenC1$values <<- eigenC1$values[1:q]
  eigenC1$vectors <<- eigenC1$vectors[,1:q]
  }))

## q first eigenpairs of C using function batchpca
system.time(replicate(10, eigenC2 <<- batchpca(C, q)))

## Check equality
all.equal(eigenC1$values, eigenC2$values)
# Reorient eigenvectors if needed
for (i in 1:q) 
{
  if (sum(eigenC1$vectors[,i] * eigenC2$vectors[,i]) < 0) 
    eigenC2$vectors[,i] <- - eigenC2$vectors[,i]
}
all.equal(eigenC1$vectors, eigenC2$vectors)
}

