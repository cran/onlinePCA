
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List sgapca_exC (const arma::colvec & ds, const arma::mat & Qs, const arma::colvec & x, const arma::colvec & gamma) {
	
	arma::mat W, R, Q = Qs;
    arma::colvec d = ds;
	arma::colvec y = trans(Q) * x;
	arma::colvec gamy = gamma % y;
	d -= gamma % d - gamy % y;
	Q += x * trans(y) * diagmat(gamma);
	qr_econ(W, R, Q);
	return Rcpp::List::create(d, W);
}



// [[Rcpp::export]]
Rcpp::List sgapca_nnC (const arma::colvec & ds, const arma::mat & Q, const arma::colvec & x, const arma::colvec & gamma) {
	
	int m = Q.n_rows, n = Q.n_cols;
    arma::colvec d = ds;
	arma::colvec y = trans(Q) * x;
	arma::colvec gamy = gamma % y;
	d -= gamma % d - gamy % y;
	arma::colvec b = y(0) * Q.col(0);
	arma::mat A(m,n);
	A.col(0) = Q.col(0) - gamy(0) * b;	
	for (int i=1; i<n; i++) {
		b += y(i-1) * Q.col(i-1) + y(i) * Q.col(i);
		A.col(i) = Q.col(i) - gamy(i) * b;
	}		
	A += x * trans(gamy);
	return Rcpp::List::create(d, A);
}
