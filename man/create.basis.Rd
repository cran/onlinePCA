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
create.basis(x, p, sp = 1e-09, degree = 3, nderiv = 2)
}

\arguments{
  \item{x}{
vector of observation times}
  \item{p}{
number of basis functions}
  \item{sp}{
smoothing parameter}
  \item{degree}{
degree of the B splines}
  \item{nderiv}{
order of the derivative to penalize for smoothing}
}

\details{
The knots of the B-spline basis are taken as regular quantiles of \code{x}. The function output is intended for use with functions 
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
	main="Sample paths of Brownian motion")

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
