sgapca <- function (d, Q, x, gamma, center, type = c("exact","nn")) 
{
    if (!missing(center)) 
    	x <- x - center
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
    k <- length(d)
	if (length(gamma) < k) 
		gamma <- rep_len(gamma,k)
	if (!is.matrix(Q))
		Q <- as.matrix(Q)
	type <- match.arg(type)
	result <- switch(type, exact = sgapca_exC(d,Q,x,gamma), 
		nn = sgapca_nnC(d,Q,x,gamma))

	ix <- order(result[[1]], decreasing = TRUE)	
    return(list(values = result[[1]][ix], vectors = result[[2]][,ix]))
}
