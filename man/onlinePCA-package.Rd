\name{onlinePCA-package}
\alias{onlinePCA-package}
\docType{package}
\title{
Online Principal Component Analysis
}
\description{
Online PCA algorithms using perturbation methods (\code{\link{perturbationRpca}}), secular equations (\code{\link{secularRpca}}), incremental PCA (\code{\link{incRpca}}), and stochastic optimization (\code{\link{ccipca}, \link{ghapca}, \link{sgapca}, \link{snlpca}}). \code{\link{impute}} handles missing data with the regression approach of Brand (2002). \code{\link{batchpca}} performs batch (offline) PCA based on covariance matrices. \code{\link{create.basis}, \link{coef2fd}, \link{fd2coef}} respectively create B-spline basis sets for functional data (FD), convert FD to basis coefficients, and convert basis coefficients back to FD. \code{\link{updateMean}} and \code{\link{updateCovariance}} update the sample mean and sample covariance.}

\details{
\tabular{ll}{
Package: \tab onlinePCA\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-6-25\cr
License: \tab GPL-3\cr
URL: \tab \url{http://cran.r-project.org/package=onlinePCA}\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
}
\author{
David Degras <ddegrasv@depaul.edu>
}
\references{
Arora et al. (2012). Stochastic Optimization for PCA and PLS.  \emph{50th Annual Conference on Communication, Control, and Computing (Allerton).}\cr
Bunch, J.R., Nielsen, C.P., and Sorensen, D.C. (1978). Rank-one modification of the symmetric eigenproblem. \emph{Numerische Mathematik.}\cr
Brand, M. (2002). Incremental singular value decomposition of uncertain data with missing values. \emph{European Conference on Computer Vision (ECCV).}\cr
Gu, M. and Eisenstat, S. C. (1994). A stable and efficient algorithm for the rank-one modification of the symmetric eigenproblem. \emph{SIAM Journal of Matrix Analysis and Applications.}\cr
Hegde et al. (2006) Perturbation-Based Eigenvector Updates for On-Line Principal Components Analysis and Canonical Correlation Analysis. \emph{Journal of VLSI Signal Processing}. \cr
Oja (1992). Principal components, Minor components, and linear neural networks. \emph{Neural Networks.}\cr
Sanger (1989). Optimal unsupervised learning in a single-layer linear feedforward neural network. \emph{Neural Networks.}\cr
Weng et al. (2003). Candid Covariance-free Incremental Principal Component Analysis. \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}. 
}
\keyword{ package }

