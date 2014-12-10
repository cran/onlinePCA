create.basis <- function(x, nknot, lambda = 1e-9, degree = 3, nderiv = 2)
{
	p <- length(x)
	knot <- quantile(x, (1:nknot)/(nknot + 1))
	delta <- c(rep.int(min(x), degree+1), knot, rep.int(max(x), degree+1))

	# B-spline design matrix
	B <- splines::splineDesign(delta, x, degree+1)
	
	# Gram matrix of B-splines
	xdiff <- diff(x, 1)
	M <- crossprod(sqrt(xdiff)*B[-1,])/2 + crossprod(sqrt(xdiff)*B[-p,])/2
	eigM <- eigen(M)
	P <- eigM$vectors
	d <- eigM$values
	sqrtM <- P %*% diag(sqrt(d)) %*% t(P)
	invsqrtM <- P %*% diag(1/sqrt(d)) %*% t(P)
	
	# Gram matrix of B-spline derivatives
	if (lambda>0)
	{
		DB <- splines::splineDesign(delta, x, degree+1, rep.int(nderiv, p))
		G <- crossprod(sqrt(xdiff)*DB[-1,])/2 + crossprod(sqrt(xdiff)*DB[-p,])/2
	}
	
	# Matrix that maps raw data to smoothed 
	# coefficients in B-spline basis
	A <- if (lambda>0) {
		crossprod(B) + (p * lambda) * G 
			} else crossprod(B)
	S <- sqrtM %*% solve(A, t(B)) 

	list(B = B, S = S, invsqrtM = invsqrtM) 
}
