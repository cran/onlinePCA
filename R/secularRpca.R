secularRpca <- function (lambda, U, x, n, f = 1/n, center, tol = 1e-10, reortho = FALSE) 
{
	if (f <= 0 || f >= 1) 
        stop("Argument 'f' must be in (0,1)")
    if (!missing(center)) 
        x <- x - center
    q <- length(lambda)
    if (q != ncol(U)) 
        stop("Arguments 'lambda' and 'U' of incompatible dimensions")

    lambda <- (1-f) * lambda
    z <- sqrt((1+f)*f) * crossprod(U,x)
    ix <- order(lambda)
    U <- U[,ix]
    lambda <- lambda[ix]
    z <- z[ix]

	# Initial deflation
    eps <- max(tol, .Machine$double.eps)
    active <- seq_len(q)
    zzero <- which(abs(z) < eps)
    if (length(zzero)) 
        active <- setdiff(active, zzero)
    ind <- which(diff(lambda[active]) < eps)
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
            U1 <- U[, g]
            z1 <- z[g]
            sigma <- sqrt(sum(z1^2))
            a <- c(sigma + z1[1], z1[-1])
            a <- a/sqrt(sum(a^2))
            H <- diag(mult) - 2 * tcrossprod(a)
            U[, g] <-U1 %*% H
            z[g] <- c(-sigma, rep(0, mult - 1))
            active <- setdiff(active, g[-1])
       }
       rm(g, U1, z1, sigma, a, H)
    }
    
    # Resolution of secular equations
    pact <- length(active)
    dact <- lambda[active]
    z2act <- z[active]^2
    bounds <- c(dact, dact[pact] + sum(z2act))
    amp <- diff(bounds)
    secular <- function(lambda) sum(z2act/{
        dact - lambda
    }) + 1
    solver <- function(i) {
        delta <- amp[i]/100
        repeat {
            lb <- bounds[i] + delta
            ub <- bounds[i + 1] - delta
            flb <- secular(lb)
            fub <- secular(ub)
            test <- flb * fub
            if (test < 0) {
                return(uniroot(secular, c(lb, ub), f.lower = flb, f.upper = fub, 
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
    lambda[active] <- roots
    
    # Computation of eigenvectors
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
    U[, active] <- U[, active] %*% eigvecs
    return(list(values = lambda[q:1], vectors = U[, q:1]))
}
