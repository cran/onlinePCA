\name{secularRpca}
\alias{secularRpca}

\title{
Recursive PCA Using Secular Equations}

\description{
The PCA is recursively updated after observation of a new vector
(rank one modification of the covariance matrix). Eigenvalues are computed as roots of a secular equation. Eigenvectors (principal components) are deduced by explicit calculation (Bunch et al., 1978) or approximated with the method of Gu and Eisenstat (1994).  
}

\usage{
secularRpca(lambda, U, x, n, f = 1/n, center, tol = 1e-10, reortho = FALSE) 
}

\arguments{
  \item{lambda}{vector of eigenvalues.} 
  \item{U}{matrix of eigenvectors (PC) stored in columns.}
 \item{x}{new data vector.} 
 \item{n}{sample size before observing \code{x}.}  
  \item{f}{forgetting factor: a number in (0,1).}
  \item{center}{optional centering vector for \code{x}.}
  \item{tol}{error tolerance for the computation of eigenvalues.}  
\item{reortho}{if FALSE, eigenvectors are explicitly computed using the method of Bunch et al. (1978). If TRUE, they are approximated with the method of Gu and Eisenstat (1994).
}
}

\details{
The method of secular equations provides accurate eigenvalues in all but pathological cases. On the other hand, the perturbation method implemented by \code{\link{perturbationRpca}} typically runs much faster but is only accurate for a large sample size \code{n}. \cr
The default eigendecomposition method is that of Bunch et al. (1978). This algorithm consists in three stages: initial deflation, nonlinear solution of secular equations, and calculation of eigenvectors.  
The calculation of eigenvectors (PCs) is accurate for the first few eigenvectors but loss of accuracy and orthogonality may occur for the next ones. In contrast the method of Gu and Eisenstat (1994) is robust against small errors in the computation of eigenvalues. It provides eigenvectors that may be less accurate than the default method but for which strict orthogonality is guaranteed. \cr
Lower values of the forgetting factor \code{f} place more weight on the current PCA while higher values place more weight on the new data. The default value \code{f = 1/n} corresponds to the usual PCA (eigenpairs of the empirical covariance matrix).\cr
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order.}
\item{vectors}{updated eigenvectors (PCs).}
}

\references{
Bunch, J.R., Nielsen, C.P., and Sorensen, D.C. (1978). Rank-one modification of the symmetric eigenproblem. \emph{Numerische Mathematik.}\cr
Gu, M. and Eisenstat, S.C. (1994). A stable and efficient algorithm for the rank-one modification of the symmetric eigenproblem. \emph{SIAM Journal of Matrix Analysis and Applications.}
}

\seealso{
\code{\link{perturbationRpca}}
}

\examples{
# Initial data set
n <- 100	
d <- 50
X <- matrix(runif(n*d),n,d)
xbar <- colMeans(X)
pca0 <- eigen(cov(X))

# New observation
newx <- runif(d)

# Recursive PCA with secular equations
xbar <- updateMean(xbar, newx, n)
secpca <- secularRpca(pca0$values, pca0$vectors, newx, n, center = xbar)

# Direct PCA
pca1 <- eigen(cov(rbind(X,newx)))

# Comparison of methods for eigenvalues
all.equal(secpca$values, pca1$values)	# TRUE

# Comparison of methods for eigenvectors
error <- sqrt(2 * abs(1 - abs(colSums(secpca$vectors * pca1$vectors ))))
plot(error, xlab = "PC", ylab = "Relative error", 
	main = "Recursive calculation of PCs\n using secular equations")
}

\keyword{multivariate}