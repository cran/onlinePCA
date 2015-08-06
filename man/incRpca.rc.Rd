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
\code{\link{incRpca}}
}
\examples{

# Data generation
n <- 300 # number of units
d <- 30000 # number of variables
n0 <- 250 # initial sample
q <- 20 # required number of PCs
X <- matrix(rnorm(n*d,sd=1/sqrt(d)),n,d) # data matrix
X <- t(apply(X,1,cumsum)) # standard Brownian motion

# Initial PCA: use duality result to transpose computations
# from dimension d * d to dimension n0 * n0 
A = scale(X[1:n0,],scale=FALSE) # center the initial data
mu0 <- attr(A,"scaled:center")
eig0 <- eigen(tcrossprod(A) / (n0-1), TRUE) # duality trick
length(eig0$values) <- q # truncation 
eig0$vectors <- eig0$vectors[,1:q]
rescale <- sqrt((n0-1)*eig0$values) # rescaling
eig0$vectors <- crossprod(A,scale(eig0$vectors, scale = rescale)) # back to dimension d 

# Incremental PCA with rotation
mu <- mu0
eig1 <- eig0
system.time({
for (i in n0:(n-1)) {
	mu <- updateMean(mu, X[i+1,], i)
	eig1 <- incRpca(eig1$values, eig1$vectors, X[i+1,], i, center = mu)
	}
})

# Incremental PCA without rotation
mu <- mu0
eig2 <- list(values = eig0$values, Ut = eig0$vectors, Us = diag(q))

system.time({
for (i in n0:(n-1)) {
	mu <- updateMean(mu, X[i+1,], i)
	eig2 <- incRpca.rc(eig2$values, eig2$Ut, eig2$Us, X[i+1,], i, center = mu)
	# rotate the PC basis and reduce its size to q every k observations
	if (i \%\% q == 0 || i == n-1) 
		{ eig2$values <- eig2$values[1:q]
		  eig2$Ut <- eig2$Ut \%*\% eig2$Us[,1:q]
		  eig2$Us <- diag(q)
	}
}
})

# Comparison
# relative differences in eigenvalues
range(eig1$values/eig2$values-1)
# (absolute) cosines of angles between eigenvectors 
abs(colSums(eig1$vectors * eig2$Ut))


}

\keyword{ multivariate }
