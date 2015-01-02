\name{perturbationRpca}
\alias{perturbationRpca}

\title{
Recursive PCA using a rank 1 perturbation method}

\description{
This function recursively updates the PCA with respect to a single new data vector, using the (fast) perturbation method of Hegde et al. (2006).}

\usage{perturbationRpca(d, Q, x, n, ff, center)}

\arguments{
  \item{d}{vector of eigenvalues}
  \item{Q}{matrix of principal components (eigenvectors of the covariance matrix) stored in columns}
  \item{x}{new data vector}
  \item{n}{current sample size (prior to observation of new data)}
  \item{ff}{forgetting factor: a number between 0 and 1, set to \code{1/n} by default}
\item{center}{centering vector (optional)}
}

\details{
The forgetting factor \code{ff} determines the balance between past and present observations in the PCA update: the closer it is to 1 (resp. to 0), the more weight is placed on current (resp. past) observations. At least one of the arguments \code{n} and \code{ff} must be specified. If \code{ff} is specified, its value overrides the argument \code{n}; otherwise, \code{ff} is set to \code{1/n} which corresponds to the assumption of a stationnary observation process. If \code{center} is not specified, then the data \code{x} is not centered prior to PCA. 
Missing values in \code{x} are imputed using the BLUP method of Brand (2002). 
}

\value{
A list with components 
  \item{values}{the updated eigenvalues.}
  \item{vectors}{the updated eigenvectors.}
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
p <- 10
X <- matrix(runif(N*p), N, p)
X <- X * rep(sqrt(12 * (1:p)), each = N)
# The eigenvalues of cov(X) are approximately equal to 1, 2, ..., p
# and the corresponding eigenvectors are approximately equal to 
# the canonical basis of R^p

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