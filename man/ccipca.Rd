\name{ccipca}
\alias{ccipca}

\title{Candid Covariance-Free Incremental PCA}

\description{Stochastic gradient ascent algorithm CCIPCA of Weng et al. (2003). 
}

\usage{
ccipca(lambda, U, x, n, q = length(lambda), l=2, center, tol = 1e-8, sort = TRUE)
}

\arguments{
\item{lambda}{vector of eigenvalues.}  
\item{U}{matrix of eigenvectors (PC) stored in columns.}
\item{x}{new data vector.}
\item{n}{sample size before observing \code{x}.} 
\item{q}{number of eigenvectors to compute.}
\item{l}{'amnesic' parameter.}
\item{center}{optional centering vector for \code{x}.}
\item{tol}{numerical tolerance.}
\item{sort}{Should the new eigenpairs be sorted?}
}

\details{
The 'amnesic' parameter \code{l} determines the weight of past observations in the PCA update. If \code{l=0}, all observations have  equal weight, which is appropriate for stationary processes. Otherwise, typical values of \code{l} range between 2 and 4. 
As \code{l} increases, more weight is placed on new observations and less on older ones. For meaningful results, the condition \code{0<=l<n} should hold.\cr
The CCIPCA algorithm iteratively updates the PCs while deflating  \code{x}. If at some point the Euclidean norm of \code{x} becomes less than \code{tol}, the algorithm stops to prevent numerical overflow.\cr
If \code{sort} is TRUE, the updated eigenpairs are sorted by decreasing eigenvalue. If FALSE, they are not sorted.      
}

\value{
A list with components  
\item{values}{updated eigenvalues.}
\item{vectors}{updated eigenvectors.}
}

\references{
Weng et al. (2003). Candid Covariance-free Incremental Principal Component Analysis. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}.
}




\examples{
## Simulation of Brownian motion
n <- 100 # number of paths
d <- 50	 # number of observation points
q <- 10	 # number of PCs to compute
x <- matrix(rnorm(n*d,sd=1/sqrt(d)), n, d)
x <- t(apply(x,1,cumsum))	

## Initial PCA
n0 <- 50
pca <- princomp(x[1:n0,])
xbar <- pca$center
pca <- list(values=pca$sdev^2, vectors=pca$loadings)

## Incremental PCA
for (i in n0:(n-1))
{	xbar <- updateMean(xbar, x[i+1,], i)
  	pca <- ccipca(pca$values, pca$vectors, x[i+1,], i, q = q, center = xbar) }

# Uncentered PCA
nx1 <- sqrt(sum(x[1,]^2))
pca <- list(values=nx1^2, vectors=as.matrix(x[1,]/nx1))
for (i in n0:(n-1))
  	pca <- ccipca(pca$values, pca$vectors, x[i+1,], i, q = q)

}