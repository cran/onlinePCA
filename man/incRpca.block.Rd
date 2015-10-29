\name{incRpca.block}
\alias{incRpca.block}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Incremental PCA with Block Update
}
\description{Sequential Karhunen-Loeve (SKL) algorithm of Levy and Lindenbaum (2000). The PCA can be updated with respect to a data matrix (not just a data vector). 
}

\usage{
incRpca.block(x, B, lambda, U, n0 = 0, f, q = length(lambda), center, byrow = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data matrix}
  \item{B}{block size}
  \item{lambda}{initial eigenvalues (optional)}
  \item{U}{initial eigenvectors/PCs (optional)}
  \item{n0}{initial sample size (optional)}
  \item{f}{vector of forgetting factors}
  \item{q}{number of requested PCs}
  \item{center}{centering vector for \code{x} (optional)}
  \item{byrow}{Are the data vectors in \code{x} stored in rows (TRUE) or columns (FALSE)?}
}
\details{
This incremental PCA algorithm utilizes QR factorization and RSVD. 
It generalizes the algorithm \code{\link{incRpca}} from vector to matrix/block updates. However, \code{\link{incRpca}} should be preferred for vector updates as it is faster.\cr
If \code{lambda} and \code{U} are specified, they are taken as the initial PCA. Otherwise, the PCA is initialized by the SVD of the first block of data in \code{x} (\code{B} data vectors). The number \code{n0} is the sample size before observing \code{x} (default value = 0). The number \code{B} is the size of blocks to be used in the PCA updates. Ideally, \code{B} should be a divisor of the number of data vectors in \code{x} (otherwise the last smaller block is discarded). If \code{U} is provided, then \code{B} and \code{q} default to the number of columns (=PC) of \code{U}.\cr 
The argument \code{f} determines the relative weight of current PCA and new data in each block update. Its length should be equal to the number of blocks in \code{x}, say \eqn{nblock}. If \code{n0} and \code{B} are provided, then \code{f} defaults to \code{B/(n0+(0:(nblock-1)*B))}, i.e., the case where all data points have equal weights. The values in \code{f} should be in (0,1), with higher values giving more weight to new data and less to the current PCA.  
}

\value{
A list with components
\item{values}{first \code{q} eigenvalues.}
\item{vectors}{first \code{q} eigenvectors/PCs.}
}

\references{
Levy, A. and Lindenbaum, M. (2000). Sequential Karhunen-Loeve basis extraction and its application to images. \emph{IEEE Transactions on Image Processing}.
}

\seealso{
\code{\link{incRpca},\link{incRpca.rc}}
}


\examples{
## Simulate Brownian Motion
n <- 100 # number of sample paths
d <- 50 # number of observation points
x <- matrix(rnorm(n*d,sd=1/sqrt(d)),n,d)
x <- t(apply(x,1,cumsum)) # dim(x) = c(100,50)
q <- 10 # number of PC to compute
B <- 20 # block size
n0 <- B # initial sample size (if relevant)

## PCA without initial values
res1 <- incRpca.block(t(x), B, q=q) # data vectors in columns
res2 <- incRpca.block(x, B, q=q, byrow=TRUE) # data vectors in rows
all.equal(res1,res2) # TRUE

## PCA with initial values 
svd0 <- svd(x[1:n0,], 0, n0) # use first block for initialization 
lambda <- svd0$d[1:n0]^2/n0 # initial eigenvalues
U <- svd0$v # initial PC
res3 <- incRpca.block(x[-(1:n0),], B, lambda, U, n0, q=q, byrow=TRUE) 
# run PCA with this initialization on rest of the data 
all.equal(res1,res3) # compare with previous PCA: TRUE

## Compare with function incRpca
res4 <- list(values=lambda, vectors=U)
for (i in (n0+1):n)
	res4 <- incRpca(res4$values, res4$vectors, x[i,], i-1, q=q)
B <- 1 # vector update
res5 <- incRpca.block(x[-(1:n0),], B, lambda, U, n0, q=q, byrow=TRUE)
ind <- which(sign(res5$vectors[1,]) != sign(res4$vectors[1,]))
res5$vectors[,ind] <- - res5$vectors[,ind] # align PCs (flip orientation as needed)
all.equal(res4,res5) # TRUE
}
