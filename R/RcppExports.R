# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ccipca_C <- function(lambda, U, x, n, q, l, tol) {
    .Call('_onlinePCA_ccipca_C', PACKAGE = 'onlinePCA', lambda, U, x, n, q, l, tol)
}

sgapca_exC <- function(Qs, x, y, gamma) {
    .Call('_onlinePCA_sgapca_exC', PACKAGE = 'onlinePCA', Qs, x, y, gamma)
}

sgapca_nnC <- function(Q, x, y, gamma) {
    .Call('_onlinePCA_sgapca_nnC', PACKAGE = 'onlinePCA', Q, x, y, gamma)
}

ghapca_C <- function(Q, x, y, gamma) {
    .Call('_onlinePCA_ghapca_C', PACKAGE = 'onlinePCA', Q, x, y, gamma)
}

