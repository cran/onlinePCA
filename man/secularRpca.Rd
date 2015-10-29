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
  \item{U}{matrix of eigenvectors (PCs) stored in columns.}
 \item{x}{new data vector.} 
 \item{n}{sample size before observing \code{x}.}  
  \item{f}{forgetting factor: a number in (0,1).}
  \item{center}{centering vector for \code{x} (optional).}
  \item{tol}{tolerance for the computation of eigenvalues.}  
\item{reortho}{if FALSE, eigenvectors are explicitly computed using the method of Bunch et al. (1978). If TRUE, they are approximated with the method of Gu and Eisenstat (1994).
}
}

\details{
The method of secular equations provides accurate eigenvalues in all but pathological cases. On the other hand, the perturbation method implemented by \code{\link{perturbationRpca}} typically runs much faster but is only accurate for a large sample size \code{n}. \cr
The default eigendecomposition method is that of Bunch et al. (1978). This algorithm consists in three stages: initial deflation, nonlinear solution of secular equations, and calculation of eigenvectors.  
The calculation of eigenvectors (PCs) is accurate for the first few eigenvectors but loss of accuracy and orthogonality may occur for the next ones. In contrast the method of Gu and Eisenstat (1994) is robust against small errors in the computation of eigenvalues. It provides eigenvectors that may be less accurate than the default method but for which strict orthogonality is guaranteed. \cr
The forgetting factor \code{f} can be interpreted as the inverse of the number of observation vectors effectively used in the PCA: the "memory" of the PCA algorithm goes back \code{1/f} observations in the past. For larger values of \code{f}, the PCA update gives more relative weight to the new data \code{x} and less to the current PCA (\code{lambda,U}). For nonstationary processes, \code{f} should be closer to 1.\cr
Only one of the arguments \code{n} and \code{f} needs being specified. If it is \code{n}, then \code{f} is set to {1/n} by default (usual PCA of sample covariance matrix where all data points have equal weight). If \code{f} is specified, its value overrides any eventual specification of \code{n}. 
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
x <- matrix(runif(n*d),n,d)
xbar <- colMeans(x)
pca0 <- eigen(cov(x))

# New observation
newx <- runif(d)

# Recursive PCA with secular equations
xbar <- updateMean(xbar, newx, n)
pca <- secularRpca(pca0$values, pca0$vectors, newx, n, center = xbar)
}
