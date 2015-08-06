\name{ghapca}
\alias{ghapca}

\title{Generalized Hebbian Algorithm for PCA}

\description{Online PCA with the GHA of Sanger (1989).}

\usage{
ghapca(lambda, U, x, gamma, q = length(lambda), center, sort = TRUE)
}

\arguments{
\item{lambda}{optional vector of eigenvalues.}  
\item{U}{matrix of eigenvectors (PC) stored in columns.}
 \item{x}{new data vector.}
 \item{gamma}{vector of gain parameters.}
 \item{q}{number of eigenvectors to compute.}
\item{center}{optional centering vector for \code{x}.}
 \item{sort}{Should the new eigenpairs be sorted?}
}

\details{
The vector \code{gamma} determines the weight placed on the new data in updating each eigenvector (the first coefficient of \code{gamma} corresponds to the first eigenvector, etc). It can be specified as a single positive number or as a vector of length \code{ncol(U)}. Larger values of \code{gamma} place more weight on \code{x} and less on \code{U}. A common choice for (the components of) \code{gamma} is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. \cr
If \code{sort} is TRUE and \code{lambda} is not missing, the updated eigenpairs are sorted by decreasing eigenvalue. Otherwise, they are not sorted.    
}

\value{
A list with components  
\item{values}{updated eigenvalues or NULL.}
\item{vectors}{updated eigenvectors.}
}

\references{
Sanger (1989). Optimal unsupervised learning in a single-layer linear feedforward neural network. \emph{Neural Networks.}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sgapca}}, \code{\link{snlpca}}
}

\examples{
## Initialization
n <- 1e4  # sample size
n0 <- 5e3 # initial sample size
d <- 10   # number of variables
q <- d # number of PC
X <- matrix(runif(n*d), n, d)
X <- X * rep(sqrt(12 * (1:d)), each = n)
# The eigenvalues of X are close to 1, 2, ..., d
# and the corresponding eigenvectors are close to 
# the canonical basis of R^d

## GHA PCA
xbar <- colMeans(X[1:n0,])
pca <- batchpca(cov(X[1:n0,]))
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, X[i,], i-1)
  pca <- ghapca(pca$values, pca$vectors, X[i,], 2/i, q, xbar)
}
pca

# Compare to batch PCA 
eigen(cov(X), TRUE)
}

\keyword{multivariate}