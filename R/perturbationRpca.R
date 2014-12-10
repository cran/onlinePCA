perturbationRpca <- function (d, Q, x, n, ff, center) 
{
    if (missing(ff)) 
        ff <- 1/n else if (ff <= 0 || ff >= 1) 
			stop("'ff' must be in (0,1)")
	k <- ncol(Q)
	p <- length(x) 
    if (length(d) != k)
    	stop("'d' and 'Q' of incompatible dimensions")
    if (nrow(Q) != p)
    	stop("'Q' and 'x' of incompatible dimensions")
    if (!missing(center)) 
    	x <- x - center
    d <- (1 - ff) * d
    na <- which(is.na(x))
    if (length(na) > 0)
    		{ if (length(na) == length(x))
    			stop("x contains only NAs")
    			 A <- Q * diag(sqrt(d))
 		if (nrow(Q) - length(na) >= q)
	 		{ x[na] <- A[na,, drop = FALSE] %*% 
 				lsfit(A[-na,, drop = FALSE], x[-na], 
 					intercept = FALSE)$coefficients
 			} else {
 			svdA <- svd(A)
 			pos <- svdA$d > sqrt(.Machine$double.eps)
 				if (!any(pos))
	 				{
	 				Ainv <- tcrossprod(svdA$v[, pos, drop = FALSE] %*% 
	 					diag(1/svdA$d[pos]), svdA$u[, pos, drop = FALSE])  
	 				x[na] <- A[na,, drop = FALSE] %*% Ainv %*% x[-na]
	 				} else x[na] <- 0
 				}	
    		}
	d2 <- d * d
    z <- sqrt((1 + ff) * ff) * crossprod(Q, x)
    z2 <- z * z
	num <- tcrossprod(z)  
	den <- matrix(d + z2, k, k, byrow = TRUE) - 
		matrix(z2 + d2, k, k)
	V <- num / den
	diag(V) <- 1
	Q <- Q %*% V
    sigma2 <- .colSums(Q * Q, p, k)
	d <- (d + z2) * sigma2
    Q <- Q * rep.int(1/sqrt(sigma2), rep.int(p,k))
   list(values = d, vectors = Q)
}
