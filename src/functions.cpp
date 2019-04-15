#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' calculate weighted correlation between columns of a matrix and a given vector
//' @export
// [[Rcpp::export]]
NumericVector matWVCorr(NumericMatrix Mat, NumericVector Vec, NumericVector Vecw) {
  arma::mat m=Rcpp::as<arma::mat>(Mat);
  arma::colvec v=Rcpp::as<arma::colvec>(Vec);
  arma::colvec vw=Rcpp::as<arma::colvec>(Vecw);
  int n=m.n_cols;
  int k=m.n_rows;
  vw/=sum(vw); // normalize weight vector
  v-=dot(v,vw);
  arma::colvec v2=v%v;
  arma::vec c(n);
  for(int i=0;i<n;i++) {
    arma::colvec ic=m.col(i);
    // shift by weighted means
    ic-=dot(ic,vw);
    double nm=dot(ic % v,vw);
    ic%=ic;
    double dn=dot(ic,vw);
    dn*=dot(v2,vw);
    c(i)=nm/sqrt(dn);
  }

  return wrap(c);
}
