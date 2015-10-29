\name{snlpca}
\alias{snlpca}

\title{Subspace Network Learning PCA}

\description{Online PCA with the SNL algorithm of Oja (1992).}

\usage{
snlpca(lambda, U, x, gamma, q = length(lambda), center, 
	type = c("exact", "nn"), sort = TRUE)
}

\arguments{
\item{lambda}{optional vector of eigenvalues.}
\item{U}{matrix of eigenvectors (PC) stored in columns.}
\item{x}{new data vector.}
\item{gamma}{vector of learning rates.}
\item{q}{number of eigenvectors to compute.}
\item{center}{optional centering vector for \code{x}.}
\item{type}{algorithm implementation: "exact" or "nn" (neural network).}
\item{sort}{Should the new eigenpairs be sorted?}

}

\details{ 
The vector \code{gamma} determines the weight placed on the new data in updating each PC. For larger values of \code{gamma}, more weight is placed on \code{x} and less on \code{U}. A common choice is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. Argument \code{gamma} can be specified as a single positive number (common to all PCs) or as a vector of length \code{q}. \cr
If \code{sort} is TRUE and \code{lambda} is not missing, the updated eigenpairs are sorted by decreasing eigenvalue. Otherwise, they are not sorted.} 

\note{The Subspace Network Learning PCA can be implemented exactly or through a neural network. The latter is less accurate but much faster. Unlike the GHA and SGA algorithms, the SNL algorithm does not consistently estimate  principal components. It provides only the linear space spanned by the PCs.
}

\value{A list with components  
\item{values}{updated eigenvalues or NULL.}
\item{vectors}{updated (rotated) eigenvectors.}
}

\references{
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}
}


\seealso{
\code{\link{ghapca}, \link{sgapca}}
}

\examples{
## Initialization
n <- 1e4  # sample size
n0 <- 5e3 # initial sample size
d <- 10   # number of variables
q <- d # number of PC to compute
x <- matrix(runif(n*d), n, d)
x <- x \%*\% diag(sqrt(12*(1:d)))
# The eigenvalues of x are close to 1, 2, ..., d
# and the corresponding eigenvectors are close to 
# the canonical basis of R^d

## SNL PCA
xbar <- colMeans(x[1:n0,])
pca <- batchpca(x[1:n0,], q, center=xbar, byrow=TRUE)
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, x[i,], i-1)
  pca <- snlpca(pca$values, pca$vectors, x[i,], 1/i, q, xbar)
}
}
