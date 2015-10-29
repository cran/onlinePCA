\name{incRpca.rc}
\alias{incRpca.rc}

\title{Incremental PCA With Reduced Complexity}

\description{
The incremental PCA is computed without rotating the updated projection space (Brand, 2002; Arora et al., 2012). Specifically, PCs are specified through a matrix of orthogonal vectors \code{Ut} that spans the PC space and a rotation matrix \code{Us} such that the PC matrix is \code{UtUs}. Given a new data vector, the PCA is updated by adding one column to \code{Ut} and recalculating the low-dimensional rotation matrix \code{Us}. This reduces complexity and helps preserving orthogonality. Eigenvalues are updated as the usual incremental PCA algorithm.}
 
\usage{
incRpca.rc(lambda, Ut, Us, x, n, f = 1/n, center, tol = 1e-07)
}

\arguments{
  \item{lambda}{vector of eigenvalues.}
  \item{Ut}{matrix of orthogonal vectors stored in columns.}
  \item{Us}{rotation matrix.}
  \item{x}{new data vector.}
  \item{n}{sample size before observing \code{x}.}
  \item{f}{forgetting factor: a number in (0,1).}
  \item{center}{optional centering vector for \code{x}.}
  \item{tol}{numerical tolerance for eigenvalues.}
}

\details{
For large datasets, this algorithm is considerably faster than its counterpart \code{\link{incRpca}}, reducing the time complexity of each update from \eqn{O(qd^2)} to \eqn{O(qd + q^3)} flops with \code{d} the length of \code{x}. A consequence of not rotating the PC basis at each update is that the dimension of the PCA decomposition increases whenever a new observation vector is not entirely contained in the PC space. To keep the number of PCs and eigenvalues from getting too large, it is necessary to multiply the matrices \eqn{U_t} and \eqn{U_s} at regular time intervals so as to recover the individual PCs and retain only the largest ones.  
}

\value{
A list with components
  \item{values}{updated eigenvalues in decreasing order.}
  \item{Ut}{updated projection space.}
  \item{Us}{updated rotation matrix.}
}
\references{
Arora et al. (2012). Stochastic Optimization for PCA and PLS.  \emph{50th Annual Conference on Communication, Control, and Computing (Allerton).}\cr
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}
}


\seealso{
\code{\link{incRpca}, \link{incRpca.block}}
}
\examples{
\dontrun{
# Data generation
n <- 400 # number of units
d <- 10000 # number of variables
n0 <- 200 # initial sample
q <- 20 # required number of PCs
x <- matrix(rnorm(n*d,sd=1/sqrt(d)),n,d) # data matrix
x <- t(apply(x,1,cumsum)) # standard Brownian motion

# Initial PCA
# Initial PCA
xbar0 <- colMeans(x[1:n0,])
pca0 <- batchpca(x0c, q, center=xbar0, byrow=TRUE)

# Incremental PCA with rotation
xbar <- xbar0
pca1 <- pca0
system.time({
for (i in n0:(n-1)) {
	xbar <- updateMean(xbar, x[i+1,], i)
	pca1 <- incRpca(pca1$values, pca1$vectors, x[i+1,], i, center=xbar)
	}
})

# Incremental PCA without rotation
xbar <- xbar0
pca2 <- list(values=pca0$values, Ut=pca0$vectors, Us=diag(q))

system.time({
for (i in n0:(n-1)) {
	xbar <- updateMean(xbar, x[i+1,], i)
	pca2 <- incRpca.rc(pca2$values, pca2$Ut, pca2$Us, x[i+1,], 
		i, center = xbar)
	# Rotate the PC basis and reduce its size to q every k observations
	if (i \%\% q == 0 || i == n-1) 
		{ pca2$values <- pca2$values[1:q]
		  pca2$Ut <- pca2$Ut \%*\% pca2$Us[,1:q]
		  pca2$Us <- diag(q)
	}
}
})

# Check that the results are identical
# Relative differences in eigenvalues
range(pca1$values/pca2$values-1)
# Cosines of angles between eigenvectors 
abs(colSums(pca1$vectors * pca2$Ut))
}
}
