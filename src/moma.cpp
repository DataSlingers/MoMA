// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"

bool DEBUG = TRUE;
using namespace arma;



// [[Rcpp::export]]
int doubleMe(const int & x) {
    return x+x;
}
<<<<<<< HEAD


// [[Rcpp::export]]
arma::vec soft_thres(arma::vec x, double l)
{
    return sign(x) % max(abs(x) - l, zeros(size(x)));
};




class Prox{ 
public:
    virtual vec prox(vec x,double l);
};

class Lasso: public Prox{
private:
public:
    vec prox(vec x,double l){
        return sign(x) % max(abs(x) - l, zeros(size(x)));
    }
};



// [[Rcpp::export]]
arma::vec test(arma::vec x,double l){
    Prox *a = new Lasso;
    return a->prox(x,l);
}

=======
>>>>>>> parent of 1ab6991... Try to build class
