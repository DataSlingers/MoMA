// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "moma.h"

#define MAX(a,b) (a)>(b)?(a):(b)
#define THRES_P(x,l) (MAX(x-l,0.0)) // shrink a positive value by `l`

/////////////////
// Section 1: Prox operators
/////////////////
inline arma::vec soft_thres(const arma::vec &x, double l){
    return arma::sign(x) % arma::max(abs(x) - l, zeros(arma::size(x)));
}

class Prox{
public:
    virtual arma::vec prox(const arma::vec &x, double l)=0;
    virtual ~Prox() = default;
};

class Lasso: public Prox{
public:
    Lasso(){
        MoMALogger::debug("A Lasso prox\n");
    }
    arma::vec prox(const arma::vec &x, double l){
        return soft_thres(x,l);
    }
};

class NNLasso: public Prox{
public:
    NNLasso(){
        MoMALogger::debug("A Non-negative Lasso prox\n");
    }
    arma::vec prox(const arma::vec &x, double l){
        return arma::max(x - l, zeros(arma::size(x)));
    }
};

class SCAD: public Prox{
private:
    double gamma; // gamma_SCAD >= 2
public:
    SCAD(double g=3.7){
        MoMALogger::debug("A SCAD prox\n");
        if(g<2) 
            MoMALogger::error("Gamma for SCAD should be larger than 2!\n");
        gamma=g;
    }
    arma::vec prox(const arma::vec &x, double l){
        int n = x.n_elem;
        arma::vec z(n);
        arma::vec absx = arma::abs(x);
        arma::vec sgn = arma::sign(x);
        // arma::vec flag = (absx >2);
        for (int i = 0; i < n; i++) // Probably need vectorization
        {
            // the implementation follows Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties
            // Jianqing Fan and Runze Li, formula(2.8)
            z(i) = absx(i) > gamma * l ? absx(i) : (absx(i) > 2 * l ? //(gamma-1)/(gamma-2) * THRES_P(absx(i),gamma*l/(gamma-1)) 
                                                    ((gamma - 1) * absx(i) - gamma * l)/ (gamma - 2)
                                                    : THRES_P(absx(i),l)
                                                    );
        }
        return z%sgn;    
    }
};


class MCP: public Prox{
private:
    double gamma; // gamma_MCP >= 1
public:
    MCP(double g=3){
        MoMALogger::debug("A MC+ prox\n");

        if(g<1) 
            MoMALogger::error("Gamma for MCP should be larger than 1!\n");
        gamma=g;
    }
    arma::vec prox(const arma::vec &x, double l){
        int n = x.n_elem;
        arma::vec z(n);
        arma::vec absx = arma::abs(x);
        arma::vec sgn = arma::sign(x);

        //// Try vectorization
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
                                    : (gamma / (gamma - 1)) * THRES_P(absx(i),l);         
        }
        return z%sgn;    
    }
};
