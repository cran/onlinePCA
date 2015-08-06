#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List ccipca_C (arma::colvec lambda, arma::mat U, 
	arma::colvec x, int n, int q, double l, double tol) {

	int i, d = x.n_elem, k = lambda.n_elem;
	if (q != k) {
		U.resize(d,q);
		lambda.resize(q);
	}
	arma::colvec v(d);  
	double f = (1.0+l)/(1.0+n);
	double nrm;

	for (i=0; i<q; i++) {  

		nrm = arma::norm(x);
		if (nrm < tol) {
			lambda.tail(q-i) = (1-f) * lambda.tail(q-i);
			break;
		}			

		if (i == n) {
			lambda[i] = nrm;			
			U.col(i) = x/nrm;
			break;
		}

		v = ((1-f) * arma::as_scalar(lambda[i])) * U.col(i) + 
			(f * arma::dot(U.col(i),x)) * x;
		nrm = arma::norm(v);
		if (nrm < tol) {
			lambda[i] = 0.0;
			break; 
		} 	
		lambda[i] = nrm; 
		U.col(i) = v/nrm;
		x -= arma::dot(U.col(i), x) * U.col(i);		

	}
	
// 	typedef	std::vector<double> stdvec;
// 	stdvec lambdav = arma::conv_to< stdvec >::from(lambda); 
// 	return Rcpp::List::create(Rcpp::Named("values") = lambdav, 
// 			Rcpp::Named("vectors") = U);
	return Rcpp::List::create(Rcpp::Named("values") = lambda, 
			Rcpp::Named("vectors") = U);
}

