fd2coef <- function(fd, basis, byrow = TRUE)
{
	if (byrow) return(tcrossprod(fd, basis$S))
	return(basis$S %*% fd)
}


