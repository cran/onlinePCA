\name{incRpca}
\alias{incRpca}

\title{Incremental PCA}

\description{This function updates the PCA using the incremental method of Arora et al. (2012).}

\usage{
incRpca(d, Q, x, n, ff = 1/n, k = ncol(Q), center)
}

\arguments{
\item{d}{eigenvalues}  
\item{Q}{principal components (PCs) stored in columns}
\item{x}{new data vector}
\item{n}{sample size prior to observing \code{x}} 
\item{ff}{forgetting factor (0 \eqn{\le} \code{ff} \eqn{\le} 1)} 
\item{k}{maximal number of PCs to compute}
\item{center}{optional centering vector for \code{x}}
}

\details{
The function provides an efficient eigendecomposition of the matrix 
\code{C = (1-ff)Qdiag(D)Q'+ \sqrt{ff(1+ff)} xx'}. 
After applying the function, the number of PCs is \code{min(k, rank([Q,x]))}.
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order}
\item{vectors}{updated principal components}
}

\references{
Arora et al. (2012). Stochastic Optimization for PCA and PLS.  \emph{50th Annual Conference on Communication, Control, and Computing (Allerton).}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
## Simulation of Brownian motion
N <- 100	 # number of paths
p <- 50	 # number of observation points
k <- 10	 # number of PCs to compute
X <- matrix(rnorm(N*p,sd=1/sqrt(p)), N, p)
X <- t(apply(X,1,cumsum))	

## Initial PCA
n0 <- 50
pca <- eigen(cov(X[1:n0,]))
xbar <- colMeans(X[1:n0,])

## Incremental PCA
for (n in n0:(N-1))
{
  xbar <- updateMean(xbar, X[n+1,], n)
  pca <- incRpca(pca$values, pca$vectors, X[n+1,], n, k = k,
		center = xbar)
}

## Batch PCA
pca1 <- eigen(cov(X))

## Error in the approximation of batch PCA (p eigenpairs) 
## by incremental PCA (k<p eigenpairs) 
par(mfrow=c(1,2))
plot(abs(1 - pca$values/pca1$values[1:k]), xlab = "PC index", log = "y",
	ylab = "Relative error", main = "Eigenvalues", type = "b")
error <- sqrt(2 * abs(1 - abs(colSums(pca$vectors * pca1$vectors[,1:k] ))))
plot(error, xlab = "PC index", ylab = "Relative error", log = "y",
	main = "Eigenvectors", type = "b")

}

\keyword{multivariate}