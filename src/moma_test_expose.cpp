# include "moma_prox.h"



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_lasso(const arma::vec &x, double l)
{
    Lasso a;
    return a.prox(x,l);
};


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_scad(const arma::vec &x, double l, double g=3.7)
{
    Scad a(g);
    return a.prox(x,l);
};


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_mcp(const arma::vec &x, double l, double g=4)
{
    Mcp a(g);
    return a.prox(x,l);
};
