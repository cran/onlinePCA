\name{snlpca}
\alias{snlpca}

\title{Subspace Network Learning PCA}

\description{This function updates the PCA using the Subspace Network Learning algorithm of Oja (1992).}

\usage{
snlpca(d, Q, x, gamma, center, type = c("exact", "nn"))
}

\arguments{
\item{d}{vector of eigenvalues (optional).}  
 \item{Q}{matrix of principal components stored in columns.}
  \item{x}{vector of new data.}
   \item{gamma}{vector of gain parameters.}
 \item{center}{vector used to center \code{x} (optional).}
  \item{type}{string specifying the algorithm implementation: "exact" or "nn" (neural network).}
}

\details{
The argument \code{d} can be omitted if the interest lies in principal components (PC's) only and not eigenvalues. \cr 
The gain vector \code{gamma} determines the weight placed on the new data in updating each PC. The first coefficient of \code{gamma} corresponds to the first PC, etc.. \code{gamma} can be specified as a single positive number (recycled by the function) or as a vector of length \code{ncol(Q)}. It  For larger values of \code{gamma}, more weight is placed on \code{x} and less on \code{Q}. A common choice is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. \cr 
The Subspace Network Learning PCA can be implemented exactly or through a neural network. The latter is less accurate but much faster.}

\value{A matrix of (rotated) principal components stored in columns.}

\references{
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}
}


\note{
Unlike the Stochastic Gradient Ascent algorithm which consistently estimates the principal components (i.e., the eigenvectors of the theoretical covariance matrix), the Subspace Network Learning algorithm only provides the linear space spanned by the PCs and not the PCs themselves. 
}

\seealso{
\code{\link{sgapca}}
}

\examples{
## Initialization
n <- 1e4  # sample size
n0 <- 5e3 # initial sample size
p <- 10   # number of variables
mat <- matrix(runif(n*p), n, p)
mat <- mat * rep(sqrt(12 * (1:p)), each = n)
# The eigenvalues of mat should be close to 1, 2, ..., p
# and the corresponding eigenvectors should be close to 
# the canonical basis of R^p

## SNL PCA
xbar <- colMeans(mat[1:n0,])
pca <- batchpca(cov(mat[1:n0,]))
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, mat[i,],  i - 1)
  pca <- snlpca(pca$values, pca$vectors, mat[i,], 1 / i, xbar)
}
pca
# Corresponding projection matrix
snlP <- tcrossprod(pca$vectors)

# Compare to batch PCA 
(batch <- eigen(cov(mat), TRUE)$vectors)
batchP <- tcrossprod(batch)

# Relative distance between the two projectors (Frobenius norm)
norm(snlP - batchP, "2") / norm(batchP, "2")
}
