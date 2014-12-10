batchpca <- function(C, k)
{
	if (missing(k))
		k <- ncol(C)
	if (k <= ncol(C) / 10) {
		res <- eigs_sym(C, k) 
	} else res <- eigen(C, TRUE)
	list(values = res$values[seq_len(k)], vectors = res$vectors[,seq_len(k)])
}
