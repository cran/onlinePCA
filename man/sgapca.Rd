\name{sgapca}
\alias{sgapca}

\title{Stochastic Gradient Ascent PCA}

\description{This function updates the PCA using the Stochastic Gradient Ascent algorithm of Oja (1992).}

\usage{
sgapca(d, Q, x, gamma, center, type = c("exact", "nn"))
}

\arguments{
\item{d}{vector of eigenvalues (optional).}  
 \item{Q}{matrix of principal components stored in columns.}
  \item{x}{new data vector.}
   \item{gamma}{vector of gain parameters.}
 \item{center}{centering vector for \code{x} (optional).}
  \item{type}{string specifying the type of implementation: "exact" or "nn" (neural network).}
}

\details{
The argument \code{d} can be omitted if the interest lies in principal components only and not  eigenvalues. Otherwise, the length of \code{d} and number of columns of \code{Q} must match.\cr 
The gain vector \code{gamma} determines the weight placed on the new data in updating each principal component. The first coefficient of \code{gamma} corresponds to the first principal component, etc.. It can be specified as a single positive number (which is recycled by the function) or as a vector of length \code{ncol(Q)}. It  For larger values of \code{gamma}, more weight is placed on \code{x} and less on \code{Q}. A common choice for (the components of) \code{gamma} is of the form \code{c/n}, with \code{n} the sample size and \code{c} a suitable positive constant. \cr 
The Stochastic Gradient Ascent PCA can be implemented exactly or through a neural network. The latter is less accurate but faster.
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order.}
\item{vectors}{updated principal components.}
}

\references{
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{snlpca}}
}

\examples{
## Initialization
n <- 1e4  # sample size
n0 <- 5e3 # initial sample size
p <- 10   # number of variables
mat <- matrix(runif(n*p), n, p)
mat <- mat * rep(sqrt(12 * (1:p)), each = n)
# The eigenvalues of mat should be close to 1, 2, ..., p
# and the corresponding eigenvectors should be close to 
# the canonical basis of R^p

## SGA PCA
xbar <- colMeans(mat[1:n0,])
pca <- batchpca(cov(mat[1:n0,]))
for (i in (n0+1):n) {
  xbar <- updateMean(xbar, mat[i,],  i - 1)
  pca <- sgapca(pca$values, pca$vectors, mat[i,], 2 / i, xbar)
}
pca

# Compare to batch PCA 
eigen(cov(mat), TRUE)
}

