\name{impute}
\alias{impute}

\title{BLUP Imputation of Missing Values}

\description{Missing values of a vector are imputed by best linear unbiased prediction (BLUP) assuming a multivariate normal distribution. 
}

\usage{
impute(lambda, U, x, center, tol = 1e-07)
}

\arguments{
\item{lambda}{vector of eigenvalues of length \code{q}.}  
\item{U}{matrix of eigenvectors (principal components) of dimension \code{p * q}.}
\item{x}{vector of observations of length \code{p} with missing entries.}
\item{center}{centering vector for \code{x}. Default is zero.}
\item{tol}{tolerance in the calculation of the pseudoinverse.}
}

\details{The vector \code{x} is assumed to arise from a multivariate normal distribution with mean vector \code{center} and covariance matrix \eqn{U diag(lambda) U^T}.}

\value{The imputed vector \code{x}.}

\references{
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}
}


\examples{
set.seed(10)
lambda <- c(1,2,5)
U <- qr.Q(qr(matrix(rnorm(30),10,3)))
x <- U \%*\% diag(sqrt(lambda)) \%*\% rnorm(3) + rnorm(10, sd =.05)
x.na <- x
x.na[c(1,3,7)] <- NA
x.imputed <- impute(lambda,U,x.na)
cbind(x,x.imputed)
}
