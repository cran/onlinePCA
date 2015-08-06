\name{perturbationRpca}
\alias{perturbationRpca}

\title{
Recursive PCA using a rank 1 perturbation method}

\description{
This function recursively updates the PCA with respect to a single new data vector, using the (fast) perturbation method of Hegde et al. (2006).}

\usage{perturbationRpca(lambda, U, x, n, f = 1/n, center, sort = TRUE)}

\arguments{
  \item{lambda}{vector of eigenvalues.}
  \item{U}{matrix of eigenvectors (PC) stored in columns.}
  \item{x}{new data vector.}
  \item{n}{sample size before observing \code{x}.}
  \item{f}{forgetting factor: a number between 0 and 1.}
\item{center}{optional centering vector for \code{x}.}
\item{sort}{Should the eigenpairs be sorted?}
}

\details{
The forgetting factor \code{f} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{f} must be specified. If \code{f} is specified, its value overrides the argument \code{n}; otherwise, \code{f} is set to \code{1/n} which corresponds to a stationnary data-generating process. Missing values in \code{x} are imputed using the BLUP method of Brand (2002). \cr
If \code{sort} is TRUE and \code{lambda} is not missing, the updated eigenpairs are sorted by decreasing eigenvalue. Otherwise, they are not sorted.  
}

\value{
A list with components 
  \item{values}{updated eigenvalues.}
  \item{vectors}{updated eigenvectors.}
}

\references{
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}\cr
Hegde et al. (2006) Perturbation-Based Eigenvector Updates for On-Line Principal Components Analysis and Canonical Correlation Analysis. \emph{Journal of VLSI Signal Processing}.}

\note{
This perturbation method is based on large sample approximations. It may be highly inaccurate for small/medium sized samples and should not be used in this case.}

\seealso{
\code{\link{secularRpca}}
}

\examples{
N <- 1e4
n0 <- 5e3
d <- 10
X <- matrix(runif(N*d), N, d)
X <- X * rep(sqrt(12 * (1:d)), each = N)
# The eigenvalues of cov(X) are approximately equal to 1, 2, ..., d
# and the corresponding eigenvectors are approximately equal to 
# the canonical basis of R^d

## Perturbation-based recursive PCA
xbar <- colMeans(X[1:n0,])
pca <- batchpca(cov(X[1:n0,]))

for (n in (n0+1):N) 
{
	x <- X[n,]
	xbar <- updateMean(xbar, x,  n)
	pca <- perturbationRpca(pca$values, pca$vectors, x, n, center = xbar)
}

## Comparison with batch PCA 
pca1 <- batchpca(cov(X))
par(mfrow=c(1,2))
plot(abs(1 - pca$values/pca1$values), xlab = "PC index", log = "y",
	ylab = "Relative error", main = "Eigenvalues", type = "b")
error <- sqrt(2 * abs(1 - abs(colSums(pca$vectors * pca1$vectors))))
plot(error, xlab = "PC index", ylab = "Relative error", log = "y",
	main = "Eigenvectors", type = "b")

}

\keyword{multivariate}