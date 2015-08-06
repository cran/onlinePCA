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
coef2fd(beta, basis, byrow = TRUE)
}

\arguments{
  \item{beta}{B-spline coefficients}
  
  \item{basis}{object created by \code{\link{create.basis}}
}
  \item{byrow}{are the coefficients of each functional observation
stored in rows (TRUE) or in columns (FALSE)?}
}

\value{
A matrix of functional data stored in the same format (row or columns) as the coefficients \code{beta}.
}

\note{In view of (online or offline) functional PCA, 
the coefficients \code{beta} are left- or right- multiplied 
by \code{M^{-1/2}} (depending on their row/column format)) 
before applying the B-spline matrix \code{B}, 
with \code{M} the Gram matrix associated to \code{B}.}



\seealso{
\code{\link{create.basis}}, \code{\link{fd2coef}}
}

\examples{
n <- 100 # number of curves
d <- 500 # number of observation points
x <- (1:d)/d # observation points
p <- 50 # number of B-spline basis functions 
par(mfrow = c(1,2))

# Simulate Brownian motion
brownian <- rnorm(n*d)/sqrt(d)
dim(brownian) <- c(n,d)
brownian <- t(apply(brownian,1,cumsum))
matplot(x, t(brownian), type = "l", ylab = "y",
	main = "Sample paths of Brownian motion")

# Create B-spline basis 
mybasis <- create.basis(x, p, 1e-4)

# Compute smooth basis coefficients 
beta <- fd2coef(brownian, mybasis)

# Recover smooth functional data
brownian.smooth <- coef2fd(beta, mybasis) 
matplot(x, t(brownian.smooth), type = "l", ylab = "y",
	main = "Smoothed sample paths")
	
# Standard PCA and Functional PCA
pca <- prcomp(brownian)
fpca <- prcomp(beta)
matplot(x, pca$rotation[,1:5], type = "l", ylab = "y", 
	main = "Standard PCA")
matplot(x, coef2fd(fpca$rotation[,1:5], mybasis, FALSE), 
	type = "l", ylab = "y", main = "Functional PCA")	
}
