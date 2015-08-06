
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::NumericMatrix sgapca_exC (const arma::mat & Qs, const arma::colvec & x,
		 const arma::colvec & y, const arma::colvec & gamma) {
	
	arma::mat W, R, Q = Qs;
	arma::colvec gamy = gamma % y;
	Q += x * arma::trans(y) * arma::diagmat(gamma);
	arma::qr_econ(W, R, Q);
	return Rcpp::wrap(W);
}



// [[Rcpp::export]]
Rcpp::NumericMatrix sgapca_nnC (const arma::mat & Q, const arma::colvec & x, 
		const arma::colvec & y, const arma::colvec & gamma) {
	
	int m = Q.n_rows, n = Q.n_cols;
 	arma::colvec gamy = gamma % y;
	arma::colvec b = y(0) * Q.col(0);
	arma::mat A(m,n);
	A.col(0) = Q.col(0) - gamy(0) * b;	
	for (int i=1; i<n; i++) {
		b += y(i-1) * Q.col(i-1) + y(i) * Q.col(i);
		A.col(i) = Q.col(i) - gamy(i) * b;
	}		
	A += x * arma::trans(gamy);
	return Rcpp::wrap(A);
}



// [[Rcpp::export]]
Rcpp::NumericMatrix ghapca_C (const arma::mat & Q, const arma::colvec & x, 
		const arma::colvec & y, const arma::colvec & gamma) {
	
	int m = Q.n_rows, n = Q.n_cols;
	arma::colvec gamy = gamma % y;
	arma::colvec b = y(0) * Q.col(0);
	arma::mat A(m,n);
	A.col(0) = Q.col(0) - gamy(0) * b;	
	for (int i=1; i<n; i++) {
		b += y(i) * Q.col(i);
		A.col(i) = Q.col(i) - gamy(i) * b;
	}		
	A += x * arma::trans(gamy);
	return Rcpp::wrap(A);
}
