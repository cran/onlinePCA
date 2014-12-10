secularRpca <- function (d, Q, x, n, ff, center, tol = 1e-10, reortho = FALSE) 
{
    if (missing(ff)) 
        ff <- 1/n
    else if (ff <= 0 || ff >= 1) 
        stop("Argument 'ff' must be in (0,1)")
    if (!missing(center)) 
        x <- x - center
    p <- length(d)
    if (p != ncol(Q)) 
        stop("Arguments 'd' and 'Q' of incompatible dimensions")
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
    z <- sqrt((1 + ff) * ff) * crossprod(Q, x)
    ix <- order(d)
    Q <- Q[, ix]
    d <- d[ix]
    z <- z[ix]
    eps <- max(tol, .Machine$double.eps)
    active <- seq_len(p)
    zzero <- which(abs(z) < eps)
    if (length(zzero)) 
        active <- setdiff(active, zzero)
    ind <- which(diff(d[active]) < eps)
    if (length(ind)) {
        ind <- sort(unique(c(ind, ind + 1L)))
        ix <- which(diff(ind) > 1)
        ub <- ind[c(ix, length(ind))]
        lb <- ind[c(1, ix + 1)]
        ngroup <- length(ub)
        group <- lapply(1:ngroup, function(i) active[lb[i]:ub[i]])
        for (i in 1:ngroup) {
            g <- group[[i]]
            mult <- length(g)
            Q1 <- Q[, g]
            z1 <- z[g]
            sigma <- sqrt(sum(z1^2))
            a <- c(sigma + z1[1], z1[-1])
            a <- a/sqrt(sum(a^2))
            H <- diag(mult) - 2 * tcrossprod(a)
            Q[, g] <-Q1 %*% H
            z[g] <- c(-sigma, rep(0, mult - 1))
            active <- setdiff(active, g[-1])
       }
       rm(g, Q1, z1, sigma, a, H)
    }
    pact <- length(active)
    dact <- d[active]
    z2act <- z[active]^2
    bounds <- c(dact, dact[pact] + sum(z2act))
    amp <- diff(bounds)
    f <- function(lambda) sum(z2act/{
        dact - lambda
    }) + 1
    solver <- function(i) {
        delta <- amp[i]/100
        repeat {
            lb <- bounds[i] + delta
            ub <- bounds[i + 1] - delta
            flb <- f(lb)
            fub <- f(ub)
            test <- flb * fub
            if (test < 0) {
                return(uniroot(f, c(lb, ub), f.lower = flb, f.upper = fub, 
                  tol = tol)$root)
            }
            else if (test == 0) {
                return(ifelse(flb == 0, lb, ub))
            }
            else if (flb > 0 && lb - bounds[i] <= tol) {
                return((bounds[i] + lb)/2)
            }
            else if (fub < 0 && bounds[i + 1] - ub <= tol) {
                return((ub + bounds[i + 1])/2)
            }
            else delta <- delta/10
        }
    }
    roots <- sapply(seq_len(pact), solver)
    d[active] <- roots
    if (reortho) {
        num <- matrix(roots - rep(dact, each = pact), pact, pact)
        den <- matrix(dact - rep(dact, each = pact), pact, pact)
        den[den == 0] <- 1
        ztilde <- sqrt(apply(num/den, 2, prod))
        eigvecs <- matrix(ztilde, pact, pact)/matrix(dact - rep(roots, 
            each = pact), pact, pact)
        norms <- sqrt(.colSums(eigvecs^2, pact, pact))
        eigvecs <- eigvecs/rep(norms, each = pact)
    }
    else {
        eigvecs <- matrix(z[active], pact, pact)/matrix(dact - 
            rep(roots, each = pact), pact, pact)
        norms <- sqrt(.colSums(eigvecs^2, pact, pact))
        eigvecs <- eigvecs/rep(norms, each = pact)
    }
    Q[, active] <- Q[, active] %*% eigvecs
    return(list(values = d[p:1], vectors = Q[, p:1]))
}
