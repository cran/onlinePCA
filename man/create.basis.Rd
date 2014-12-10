\name{create.basis}
\alias{create.basis}

\title{
Create a smooth B-spline basis
}

\description{
This function creates a smooth B-spline basis and provides tools to
find the coefficients of functional data in the basis  
and to recover functional data from basis coefficients.}

\usage{
create.basis(x, nknot, lambda = 1e-09, degree = 3, nderiv = 2)
}

\arguments{
  \item{x}{
vector of observation times}
  \item{nknot}{
number of interior knots}
  \item{lambda}{
smoothing parameter}
  \item{degree}{
degree of the B splines}
  \item{nderiv}{
order of the derivative to penalize for smoothing}
}

\details{
The knots of the B-spline basis are taken as regular quantiles of \code{x}. The total number of knots is \code{nknot + degree + 1}. The function output is typically not meant to be used directly but rather through the functions 
\code{\link{coef2fd}} and \code{\link{fd2coef}}. 
}

\value{
A list with fields
\item{B}{matrix of B-splines evaluated at \code{x}
(each column represents a basis function)}
\item{S}{matrix that maps functional data
to their (smoothed) coefficients of their projection
in the basis set. For the purpose of PCA, the 
coefficients are premultiplied by \code{M^{1/2}}, 
where \code{M} is the Gram matrix associated with \code{B}}
\item{invsqrtM}{matrix \code{M^{-1/2}} used to recover functional
# data from their coefficients in the basis set}
}

\seealso{
\code{\link{coef2fd}}, \code{\link{fd2coef}}
}

\examples{
n <- 100	# number of curves
p <- 500	# number of observation points
x <- (1:p)/p	# observation grid
nknot <- 50  # number of interior knots 
par(mfrow = c(1,2))

# Simulate Brownian motion
brownian <- rnorm(n*p)/sqrt(p)
dim(brownian) <- c(n,p)
brownian <- t(apply(brownian,1,cumsum))
matplot(x, t(brownian), type = "l", ylab = "y",
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
