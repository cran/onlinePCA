\name{perturbationRpca}
\alias{perturbationRpca}

\title{
Recursive PCA using a rank 1 perturbation method}

\description{
This function recursively updates the PCA with respect to a single new data vector, using the (fast) perturbation method of Hegde et al. (2006).}

\usage{perturbationRpca(lambda, U, x, n, f = 1/n, center, sort = TRUE)}

\arguments{
  \item{lambda}{vector of eigenvalues.}
  \item{U}{matrix of eigenvectors (PC) stored in columns.}
  \item{x}{new data vector.}
  \item{n}{sample size before observing \code{x}.}
  \item{f}{forgetting factor: a number between 0 and 1.}
\item{center}{optional centering vector for \code{x}.}
\item{sort}{Should the eigenpairs be sorted?}
}

\details{
The forgetting factor \code{f} can be interpreted as the inverse of the number of observation vectors effectively used in the PCA: the "memory" of the PCA algorithm goes back \code{1/f} observations in the past. For larger values of \code{f}, the PCA update gives more relative weight to the new data \code{x} and less to the current PCA (\code{lambda,U}). For nonstationary processes, \code{f} should be closer to 1.\cr
Only one of the arguments \code{n} and \code{f} needs being specified. If it is \code{n}, then \code{f} is set to {1/n} by default (usual PCA of sample covariance matrix where all data points have equal weight). If \code{f} is specified, its value overrides any eventual specification of \code{n}. 
 \cr
If \code{sort} is TRUE, the updated eigenpairs are sorted by decreasing eigenvalue. Otherwise, they are not sorted.  
}

\value{
A list with components 
  \item{values}{updated eigenvalues.}
  \item{vectors}{updated eigenvectors.}
}

\references{
Hegde et al. (2006) Perturbation-Based Eigenvector Updates for On-Line Principal Components Analysis and Canonical Correlation Analysis. \emph{Journal of VLSI Signal Processing}.}

\note{
This perturbation method is based on large sample approximations. It tends to be highly inaccurate for small/medium sized samples and should not be used in this case.}

\seealso{
\code{\link{secularRpca}}
}

\examples{
n <- 1e3
n0 <- 5e2
d <- 10
x <- matrix(runif(n*d), n, d)
 x <- x \%*\% diag(sqrt(12*(1:d)))
# The eigenvalues of cov(x) are approximately equal to 1, 2, ..., d
# and the corresponding eigenvectors are approximately equal to 
# the canonical basis of R^d

## Perturbation-based recursive PCA
# Initialization: use factor 1/n0 (princomp) rather 
# than factor 1/(n0-1) (prcomp) in calculations
pca <- princomp(x[1:n0,], center=FALSE)
xbar <- pca$center
pca <- list(values=pca$sdev^2, vectors=pca$loadings) 

for (i in (n0+1):n) {
	xbar <- updateMean(xbar, x[i,], i-1)
	pca <- perturbationRpca(pca$values, pca$vectors, x[i,], 
		i-1, center=xbar) }
}
