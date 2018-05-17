// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"

using namespace arma;


class Prox{
public:
    virtual arma::vec thres(const arma::vec &x, double l)=0;
};

class Lasso: public Prox{
public:
    arma::vec thres(const arma::vec &x, double l){
        Rcpp::Rcout<<("In Lasso.prox\n");
        return sign(x) % max(abs(x) - l, zeros(size(x)));
    }
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec test(arma::vec x, double l)
{
    Lasso a;
    
    Rcpp::Rcout<< "In Test\n"<< std::endl;
    return a.thres(x,l);;
};
