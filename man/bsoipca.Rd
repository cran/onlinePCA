\name{bsoipca}
\alias{bsoipca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Block Stochastic Orthononal Iteration (BSOI)}
\description{
The online PCA algorithm of Mitliagkas et al. (2013) is a block-wise stochastic variant of the classical power-method.}

\usage{
bsoipca(x, q, U, B, center, byrow = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data matrix.}
  \item{q}{number of PC to compute.}
  \item{U}{matrix of initial PCs in columns (optional).}
  \item{B}{size of block updates (optional).}
  \item{center}{centering vector (optional).}
  \item{byrow}{are the data vectors in \code{x} stored in rows (TRUE) or columns (FALSE)?
}
}
\details{
The default value of \code{B} is \eqn{floor(n/nblock)} with \eqn{n} the number of data vectors in \code{x}, \eqn{d} the number of variables, and \eqn{nblock=ceiling(log(d))} the number of blocks.\cr
If \code{U} is specified, \code{q} defaults to \code{ncol(U)}; otherwise the initial PCs are computed from the first block of data and \code{q} must be specified explicitly. \cr
Although the algorithm does not give eigenvalues, they can easily be estimated by computing the variance of the data along the PCs.  
}
\value{
A matrix with the \code{q} first eigenvectors/PCs in columns.
}

\references{
Mitliagkas et al. (2013). Memory limited, streaming PCA. \emph{Advances in Neural Information Processing Systems}.}


\examples{
## Simulate Brownian Motion
n <- 100 # number of sample paths
d <- 50 # number of observation points
x <- matrix(rnorm(n*d,sd=1/sqrt(d)),n,d)
x <- t(apply(x,1,cumsum)) # dim(x) = c(100,50)

q <- 10 # number of PC to compute
B <- 20 # block size

## BSOI PCA 
U <- bsoipca(x, q, B=B, byrow=TRUE) # PCs
lambda <- apply(x \%*\% U, 2, var) # eigenvalues 
}
