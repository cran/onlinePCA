\name{onlinePCA-package}
\alias{onlinePCA-package}
\docType{package}
\title{
Online Principal Component Analysis
}
\description{
Online PCA using perturbation methods (\code{\link{perturbationRpca}, \link{secularRpca}}), incremental PCA (\code{\link{incRpca}}), and stochastic gradient approaches (\code{\link{sgapca}, \link{snlpca}}). Missing data are handled by the regression approach of Brand (2002). \code{\link{batchpca}} computes batch (offline) PCA based on covariance matrices. \code{\link{create.basis}, \link{coef2fd}, \link{fd2coef}} respectively create B-spline basis functions for functional data (FD), compute the basis coefficients of FD, and convert basis coefficients to FD. \code{\link{updateMean}} and \code{\link{updateCovariance}} update the sample mean and sample covariance.}

\details{
\tabular{ll}{
Package: \tab onlinePCA\cr
Type: \tab Package\cr
Version: \tab 1.0-1\cr
Date: \tab 2014-12-10\cr
License: \tab GPL-3\cr
URL: \tab \url{http://cran.r-project.org/package=onlinePCA}\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
David Degras <ddegrasv@depaul.edu>
}
\references{
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}\cr
Gu, M. and Eisenstat, S. C. (1994). A stable and efficient algorithm for the rank-one modification of the symmetric eigenproblem. \emph{SIAM Journal of Matrix Analysis and Applications.}\cr
Hegde et al. (2006) Perturbation-Based Eigenvector Updates for On-Line Principal Components Analysis and Canonical Correlation Analysis. \emph{Journal of VLSI Signal Processing}. \cr
Li, W., Yue, H. H., Valles-Cervantes, S. and Qin, S. J. (2000). Recursive PCA for adaptive process monitoring. \emph{Journal of Process Control}.\cr
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.} 
}
\keyword{ package }

