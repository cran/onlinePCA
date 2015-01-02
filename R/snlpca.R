snlpca <- function (d, Q, x, gamma, center, type = c("exact", "nn")) 
{
	if (!missing(center)) 
		x <- x - center
	type <- match.arg(type)	
    y <- as.numeric(crossprod(Q, x))
    gamy <- gamma * y
	d <- if (missing(d)) NULL else (1 - gamma) * d + gamy * y
    na <- which(is.na(x))
    if (length(na) > 0)
    		{ if (length(na) == length(x))
    			stop("x contains only NAs")
    			 A <- Q %*% diag(sqrt(d))
 		if (nrow(Q) - length(na) >= q)
	 		{ ginvAx.nona <- suppressWarnings(lsfit(A[-na,, drop = FALSE], 
	 			x[-na], intercept = FALSE)$coefficients)
	 			x[na] <- A[na,, drop = FALSE] %*% ginvAx.nona
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
	
	if (type == "exact") {
	    Q <- Q + tcrossprod(x, gamy)
	    eig <- eigen(crossprod(Q), TRUE)
	    nonzero <- which(eig$values > sqrt(.Machine$double.eps))
	    iS <- eig$vectors[, nonzero] %*% 
	    	(t(eig$vectors[, nonzero])/sqrt(eig$values[nonzero]))
	    Q <- Q %*% iS
	if (length(nonzero) < ncol(Q)) 
        warning(paste("Matrix 'Q' is not full rank. Returning", 
            length(nonzero), "PC."))
	} else if (type == "nn") {
		Q <- Q + tcrossprod(x - Q %*% y, gamy)
	}

	return(list(values = d, vectors = Q))
}
