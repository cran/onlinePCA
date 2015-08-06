batchpca <- function(C, q)
{
	if (missing(q))
		q <- ncol(C)
	if (q <= ncol(C) / 10) {
		res <- eigs_sym(C, q) 
	} else res <- eigen(C, TRUE)
	list(values = res$values[seq_len(q)], vectors = res$vectors[,seq_len(q)])
}
