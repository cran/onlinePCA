incRpca <- function (d, Q, x, n, ff = 1/n, k = ncol(Q), center) 
{
    q <- length(d)
    if (ncol(Q) != q) 
        stop("length(d) != ncol(Q)")
    if (nrow(Q) != length(x)) 
        stop("length(x) != nrow(Q)")
    if (!missing(center)) 
        x <- x - center
    	d <- (1 - ff) * d
    	x <- sqrt(ff * (1 + ff)) * x
    na <- which(is.na(x))
    if (length(na) > 0)
    		{ if (length(na) == length(x))
    			stop("x contains only NAs")
    			 A <- Q %*% diag(sqrt(d))
 		if (nrow(Q) - length(na) >= q)
	 		{ ginvAx.nona <- suppressWarnings(lsfit(A[-na,, drop = FALSE], x[-na], intercept = FALSE)$coefficients)
	 			x[na] <- A[na,, drop = FALSE] %*% ginvAx.nona
 			} else {
 			svdA <- svd(A)
 			pos <- svdA$d > sqrt(.Machine$double.eps)
 				if (!any(pos))
	 				{
	 				Ainv <- tcrossprod(svdA$v[, pos, drop = FALSE] %*% diag(1/svdA$d[pos]), svdA$u[, pos, drop = FALSE])  
	 				x[na] <- A[na,, drop = FALSE] %*% Ainv %*% x[-na]
	 				} else x[na] <- 0
 				}	
    		}
    xhat <- crossprod(Q, x)
    xorth <- x - Q %*% xhat
    xorth_ <- sqrt(sum(xorth^2))
    M <- matrix(nrow = q + 1L, ncol = q + 1L)
    M[1:q, 1:q] <- diag(d) + tcrossprod(xhat)
    M[1:q, q + 1L] <- xorth_ * xhat
    M[q + 1L, 1:q] <- xorth_ * xhat
    M[q + 1L, q + 1L] <- xorth_^2
    eigM <- eigen(M, TRUE)
    k <- min(k, q + 1L - (xorth_ < sqrt(.Machine$double.eps)))
    Q <- cbind(Q, xorth/xorth_) %*% eigM$vectors[, 1:k]
    return(list(values = eigM$values[1:k], vectors = Q))
}

