\name{fd2coef}
\alias{fd2coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the coefficients of functional data in a B-spline basis
}

\description{
This function computes the coefficients of functional data in a B-spline basis.
}

\usage{
fd2coef(fd, basis, byrow = TRUE)
}

\arguments{
  \item{fd}{matrix of functional data.
}
  \item{basis}{
object created by \code{\link{create.basis}}}.
  \item{byrow}{
are the functional data stored in rows (TRUE) or in columns (FALSE)?}
}

\value{
A matrix of B-spline coefficients stored in the same format (row or columns) as functional data.
}

\note{In view of (online or offline) functional PCA, 
the coefficients are smoothed and premultiplied by \code{M^{1/2}}, 
with \code{M} the Gram matrix associated to the B-spline matrix.}

\seealso{
\code{\link{create.basis}}, \code{\link{coef2fd}}
}

\examples{
n <- 100	# number of curves
p <- 500	# number of observation points
x <- (1:p)/p	# observation grid
nknot <- 50  # number of interior knots 
par(mfrow=c(1,2))

# Simulate Brownian motion
brownian <- rnorm(n*p)/sqrt(p)
dim(brownian) <- c(n,p)
brownian <- t(apply(brownian,1,cumsum))
matplot(x, t(brownian), type="l", ylab="y",
	main="Sample paths of Brownian motion")

# Create B-spline basis 
mybasis <- create.basis(x, nknot, 1e-4)

# Compute smooth basis coefficients 
theta <- fd2coef(brownian, mybasis)

# Recover smooth functional data
brownian.smooth <- coef2fd(theta, mybasis) 
matplot(x, t(brownian.smooth), type = "l", ylab = "y",
	main = "Smoothed sample paths")
	
# Standard PCA and Functional PCA
pca <- prcomp(brownian)
fpca <- prcomp(theta)
matplot(x, pca$rotation[,1:5], type = "l", ylab = "y", 
	main = "Standard PCA")
matplot(x, coef2fd(fpca$rotation[,1:5], mybasis, FALSE), 
	type = "l", ylab = "y", main = "Functional PCA")

	
}
