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
fd2coef(x, basis, byrow = TRUE)
}

\arguments{
  \item{x}{matrix of functional data.}
  \item{basis}{object created by \code{\link{create.basis}}}.
  \item{byrow}{are the functional data stored in rows (TRUE) or in columns (FALSE)?}
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
n <- 100 # number of curves
d <- 500 # number of observation points
grid <- (1:d)/d # observation points
p <- 50 # number of B-spline basis functions 

# Simulate Brownian motion
x <- matrix(rnorm(n*d,sd=1/sqrt(d)),n,d)
x <- t(apply(x,1,cumsum))

# Create B-spline basis 
mybasis <- create.basis(grid, p, 1e-4)

# Compute smooth basis coefficients 
beta <- fd2coef(x, mybasis)

# Recover smooth functional data
x.smooth <- coef2fd(beta, mybasis) 
	
# Standard PCA and Functional PCA
pca <- prcomp(x)
fpca <- prcomp(beta)
}
