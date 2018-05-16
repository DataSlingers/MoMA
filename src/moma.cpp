// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"

using namespace arma;


class Prox{
public:
    virtual arma::vec thres(arma::vec x, double l);
};

class Lasso: public Prox{
public:
    arma::vec thres(arma::vec x, double l){
        printf("In Lasso.prox");
        return sign(x) % max(abs(x) - l, zeros(size(x)));
    }
};

// [[Rcpp::export]]
arma::vec test(arma::vec x, double l)
{
    Lasso a;
// not working
    return x;
}
