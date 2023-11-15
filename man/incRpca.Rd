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
The forgetting factor \code{f} can be interpreted as the inverse of the number of observation vectors effectively used in the PCA: the "memory" of the PCA algorithm goes back \code{1/f} observations in the past. For larger values of \code{f}, the PCA update gives more relative weight to the new data \code{x} and less to the current PCA (\code{lambda,U}). For nonstationary processes, \code{f} should be closer to 1.\cr
Only one of the arguments \code{n} and \code{f} needs being specified. If it is \code{n}, then \code{f} is set to \code{1/n} by default (usual PCA of sample covariance matrix where all data points have equal weight). If \code{f} is specified, its value overrides any eventual specification of \code{n}. 
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
## Simulate Brownian motion
n <- 100 # number of sample paths
d <- 50	 # number of observation points
q <- 10	 # number of PCs to compute
n0 <- 50 # number of sample paths used for initialization 
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))	
dim(x) # (100,50)


## Incremental PCA (IPCA, centered)
pca <- prcomp(x[1:n0,]) # initialization
xbar <- pca$center
pca <- list(values=pca$sdev[1:q]^2, vectors=pca$rotation[,1:q])
for (i in (n0+1):n)
{
  xbar <- updateMean(xbar, x[i,], i-1)
  pca <- incRpca(pca$values, pca$vectors, x[i,], i-1, q = q,
		center = xbar)
}

## Incremental PCA (IPCA, uncentered)
pca <- prcomp(x[1:n0,],center=FALSE) # initialization
pca <- list(values = pca$sdev[1:q]^2, vectors = pca$rotation[,1:q])
for (i in (n0+1):n)
  pca <- incRpca(pca$values, pca$vectors, x[i,], i-1, q = q)
}
