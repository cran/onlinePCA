\name{incRpca}
\alias{incRpca}

\title{Incremental PCA}

\description{
Online PCA using the incremental SVD method of Brand (2002) and Arora et al. (2012).}

\usage{
incRpca(lambda, U, x, n, f = 1/n, q = length(lambda), center, tol = 1e-7)
}

\arguments{
\item{lambda}{vector of eigenvalues.}  
\item{U}{matrix of eigenvectors (principal components) stored in columns.}
\item{x}{new data vector.}
\item{n}{sample size before observing \code{x}.} 
\item{f}{forgetting factor: a number in (0,1).} 
\item{q}{number of eigenvectors to compute.}
\item{center}{optional centering vector for \code{x}.}
\item{tol}{numerical tolerance.}
}

\details{
If the Euclidean distance between \code{x} and \code{U} is more than \code{tol}, the number of eigenpairs increases to \code{length(lambda)+1} before eventual truncation at order \code{q}. Otherwise, the eigenvectors remain unchanged and only the eigenvalues are updated.
\cr  
The forgetting factor \code{f} can be interpreted as the inverse of the number of observation vectors efectively used in the PCA: the "memory" of the PCA algorithm goes back \code{1/f} observations in the past. The default value \code{f=1/n} corresponds to a stationary process (incremental PCA of the empirical covariance matrix). For larger values of \code{f}, the PCA update places more weight on the new observation \code{x} (and less on the current PCA) in the update. For nonstationary processes, \code{f} should be closer to 1.  
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order.}
\item{vectors}{updated eigenvectors.}
}

\references{
Arora et al. (2012). Stochastic Optimization for PCA and PLS.  \emph{50th Annual Conference on Communication, Control, and Computing (Allerton).}\cr
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}
}




\examples{
## Simulation of Brownian motion
N <- 100	 # number of paths
d <- 50	 # number of observation points
q <- 10	 # number of PCs to compute
X <- matrix(rnorm(N*d,sd=1/sqrt(d)), N, d)
X <- t(apply(X,1,cumsum))	

## Initial PCA
n0 <- 50
pca <- eigen(cov(X[1:n0,]))
xbar <- colMeans(X[1:n0,])

## Incremental PCA
for (n in n0:(N-1))
{
  xbar <- updateMean(xbar, X[n+1,], n)
  pca <- incRpca(pca$values, pca$vectors, X[n+1,], n, q = q,
		center = xbar)
}

## Batch PCA
pca1 <- eigen(cov(X))

## Error in the approximation of batch PCA (p eigenpairs) 
## by incremental PCA (q<d eigenpairs) 
par(mfrow=c(1,2))
plot(abs(1 - pca$values/pca1$values[1:q]), xlab = "PC index", log = "y",
	ylab = "Relative error", main = "Eigenvalues", type = "b")
error <- sqrt(2 * abs(1 - abs(colSums(pca$vectors * pca1$vectors[,1:q] ))))
plot(error, xlab = "PC index", ylab = "Relative error", log = "y",
	main = "Eigenvectors", type = "b")

}

\keyword{multivariate}