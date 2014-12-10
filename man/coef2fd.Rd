\name{coef2fd}
\alias{coef2fd}

\title{
Recover functional data from their B-spline coefficients
}
\description{
This function computes functional data 
from their coefficients in a B-spline basis.
}

\usage{
coef2fd(theta, basis, byrow = TRUE)
}

\arguments{
  \item{theta}{B-spline coefficients}
  
  \item{basis}{object created by \code{\link{create.basis}}
}
  \item{byrow}{are the coefficients of each functional observation
stored in rows (TRUE) or in columns (FALSE)?}
}

\value{
A matrix of functional data stored in the same format (row or columns) as  
the coefficients \code{theta}.
}

\note{In view of (online or offline) functional PCA, 
the coefficients \code{theta} are left- or right- multiplied 
by \code{M^{-1/2}} (depending on their row/column format)) 
before applying the B-spline matrix \code{B}, 
with \code{M} the Gram matrix associated to \code{B}.}



\seealso{
\code{\link{create.basis}}, \code{\link{fd2coef}}
}

\examples{
n <- 100 # number of curves
p <- 500	 # number of observation points
x <- (1:p)/p # observation grid
nknot <- 50 # number of interior knots 
par(mfrow = c(1,2))

# Simulate Brownian motion
brownian <- rnorm(n*p)/sqrt(p)
dim(brownian) <- c(n,p)
brownian <- t(apply(brownian,1,cumsum))
matplot(x, t(brownian), type = "l", ylab = "y",
	main = "Sample paths of Brownian motion")

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
