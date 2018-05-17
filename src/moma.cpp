// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include "moma.h"
using namespace arma;
#define MAX(a,b) (a)>(b)?(a):(b)


 class Prox{
public:
    virtual arma::vec prox(const arma::vec &x, double l)=0;
};

class Lasso: public Prox{
public:
    arma::vec prox(const arma::vec &x, double l){
        Rcpp::Rcout<<("In Lasso.prox\n");
        return sign(x) % arma::max(abs(x) - l, zeros(size(x)));
    }
};

class Scad: public Prox{
private:
    double gamma;
    
public:
    Scad(double g=3.7){gamma=g;}
    arma::vec prox(const arma::vec &x, double l){
        int n = x.n_elem;
        arma::vec z(n);
        arma::vec abs = arma::abs(x);
        arma::vec sgn = sign(x);
        for (int i = 0; i < n; i++) // Probably need vectorization
        {
            z(i) = abs(i) > gamma * l ? x(i)
                                    : (abs(i) > 2 * l ? ((gamma - 1) * x(i) - sgn(i) * gamma * l) / (gamma - 2)
                                                        : sgn(i) * MAX(double(abs(i) - l),0.0));
        }
        return z;    
    }
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec test(arma::vec x, double l)
{
    Lasso a;
    Rcpp::Rcout<< "In Test\n"<< std::endl;
    return a.prox(x,l);;
};
