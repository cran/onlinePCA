\name{secularRpca}
\alias{secularRpca}

\title{
Recursive PCA using secular equations}

\description{
This function recursively updates a PCA with respect to a new data vector. The eigenvalues are computed as the roots of a secular equations. The principal components are either deduced by direct calculation or computed with the approximation method of Gu and Eisenstat (1994).  
}

\usage{
secularRpca(d, Q, x, n, ff, center, tol = 1e-10, reortho = FALSE) 
}

\arguments{
  \item{d}{vector of eigenvalues} 
  \item{Q}{matrix of principal components (PCs) stored in columns}
 \item{x}{new data vector} 
 \item{n}{sample size before observing \code{x}}  
  \item{ff}{forgetting factor; a real number between 0 and 1 set to \code{1/n} by default.}
  \item{center}{optional centering vector for \code{x}}
  \item{tol}{error tolerance for the computation of eigenvalues}  
\item{reortho}{logical; if FALSE (default) 
	. If TRUE, the approximation method of Gu and Eisenstat (1994) is used to compute the eigenvectors. 
}
}

\details{
The method of secular equations provides accurate eigenvalues in all but pathological cases. On the other hand, the perturbation method implemented by \code{\link{perturbationRpca}} typically runs much faster but is only accurate for a large sample size \code{n}. \cr
The default method for computing PCs/eigenvalues (\code{reortho = FALSE}) is that of Li et al. (2000). It is accurate for the first few PCs but loss of accuracy and orthogonality may occur for the next PCs. In contrast the method of Gu and Eisenstat (1994) is robust against small errors in the computation of eigenvalues. It provides PCs that may be less accurate than the default method  but for which strict orthogonality is guaranteed.\cr
Lower values of the forgetting factor \code{ff} place more weight on the current PCA while higher values place more weight on the new data. The default value \code{ff = 1/n} corresponds to the usual PCA (eigenpairs of the empirical covariance matrix).\cr
Missing values in \code{x} are imputed using the BLUP method of Brand (2002). 
}

\value{
A list with components  
\item{values}{updated eigenvalues in decreasing order}
\item{vectors}{updated principal components}
}

\references{
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}\cr
Gu, M. and Eisenstat, S. C. (1994). A stable and efficient algorithm for the rank-one modification of the symmetric eigenproblem. \emph{SIAM Journal of Matrix Analysis and Applications.}\cr
Li, W., Yue, H. H., Valles-Cervantes, S. and Qin, S. J. (2000). Recursive PCA for adaptive process monitoring. \emph{Journal of Process Control}.
}

\seealso{
\code{\link{perturbationRpca}}
}

\examples{
# Initial data set
n <- 100	
p <- 50
X <- matrix(runif(n*p),n,p)
xbar <- colMeans(X)
pca0 <- eigen(cov(X))

# New observation
newx <- runif(p)

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