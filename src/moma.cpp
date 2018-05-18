// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
#include "moma_logging.h"
#include "moma.h"

using namespace Rcpp;
using namespace arma;
#define MAX(a,b) (a)>(b)?(a):(b)

inline arma::vec soft_thres(const arma::vec &x, double l){
    return sign(x) % arma::max(abs(x) - l, zeros(size(x)));
}
inline double soft_thres_e(double x, double l){
    double sgn;
    if(x>0) sgn = 1;
    else sgn = -1;
    return sgn * MAX(sgn*x - l,0);
}
class Prox{
public:
    virtual arma::vec prox(const arma::vec &x, double l)=0;
};

class Lasso: public Prox{
public:
    arma::vec prox(const arma::vec &x, double l){
        return soft_thres(x,l);
    }
};

class Scad: public Prox{
private:
    double gamma; // >= 2
    
public:
    Scad(double g=3.7){
        if(g<2) 
        Rcpp::stop("Gamma for MCP should be larger than 2!\n");
        gamma=g;
    }
    arma::vec prox(const arma::vec &x, double l){
        int n = x.n_elem;
        arma::vec z(n);
        arma::vec absx = arma::abs(x);
       
        arma::vec sgn = sign(x);
        for (int i = 0; i < n; i++) // Probably need vectorization
        {
            // the implementation follows Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties
            // Jianqing Fan and Runze Li, formula(2.8)

            z(i) = absx(i) > gamma * l ? absx(i) : (
                                                    absx(i) > 2 * l ? ((gamma - 1) * absx(i) - gamma * l)/ (gamma - 2)
                                                   // absx(i) > 2 * l ? (gamma-1)/(gamma-2)*sgn(i) * MAX(double(absx(i) - gamma*l/(gamma-1)),0.0)
                                                    : MAX(double(absx(i) - l),0.0)
                                                );
        }
        return z%sgn;    
    }
};


class Mcp: public Prox{
private:
    double gamma; // >= 1
public:
    Mcp(double g=4){
        if(g<1) stop("Gamma for MCP should be larger than 1!\n");
        gamma=g;
    }
    arma::vec prox(const arma::vec &x, double l){
        int n = x.n_elem;
        arma::vec z(n);
        arma::vec absx = arma::abs(x);
        arma::vec sgn = arma::sign(x);

        // arma::vec thr = sgn % arma::max(absx - l, zeros(size(x)));
        // arma::vec flag = ones<vec>(n) * gamma*l;
        // arma::vec large = x>flag;
        // arma::vec small = ones(gamma*l)-large;
        for (int i = 0; i < n; i++) // Probably need vectorization
        {
            // implementation follows lecture notes of Patrick Breheny
            // http://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-29.pdf
            // slide 19
            z(i) = absx(i) > gamma * l ? absx(i)
                                    : (gamma/(gamma-1)) * double(MAX(double(absx(i) - l),0.0));
            // if(absx(i) <= gamma*l){
            //     Rcout << i << "\t";
            //     Rcout << sgn(i)*(gamma/(gamma-1)) <<"\t";
            //     Rcout << (MAX(double(absx(i) - l),0.0));
            //     Rcout << endl;
            // }
           
        }
        return z%sgn;    
    }
};

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_lasso(arma::vec x, double l)
{
    Lasso a;
    return a.prox(x,l);
};


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_scad(arma::vec x, double l, double g=3.7)
{
    Scad a(g);
    return a.prox(x,l);
};


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec prox_mcp(arma::vec x, double l, double g=4)
{
  
    Mcp a(g);
    return a.prox(x,l);
};
