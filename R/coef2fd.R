coef2fd <- function(theta, basis, byrow = TRUE)
{
	if (byrow) return(theta %*% basis$invsqrtM %*% t(basis$B))
	return(basis$B %*% basis$invsqrtM %*% theta)	
}
