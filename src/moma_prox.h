// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef MOMA_PROX
#define MOMA_PROX 1

#include "moma_base.h"
#include "moma_logging.h"

#define MAX(a,b) (a)>(b)?(a):(b)
#define THRES_P(x,l) (MAX(x-l,0.0)) // shrink a positive value by `l`

inline arma::vec soft_thres(const arma::vec &x, double l){
    return arma::sign(x) % arma::max(abs(x) - l, zeros(arma::size(x)));
}

class Prox{
public:
    virtual arma::vec prox(const arma::vec &x, double l) = 0;
    virtual ~Prox() = default;
};

class Lasso: public Prox{
public:
    Lasso(){
        MoMALogger::debug("Initializing Lasso proximal operator object");
    }
    arma::vec prox(const arma::vec &x, double l){
        return soft_thres(x,l);
    }
};

class SCAD: public Prox{
private:
    double gamma; // gamma_SCAD >= 2
public:
    SCAD(double g = 3.7){
        MoMALogger::debug("Initializing SCAD proximal operator object");

        if(g < 2){
            MoMALogger::error("Non-convexity parameter (gamma) must be larger than or equal to 2 for SCAD");
        }

        gamma = g;
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
        return z % sgn;
    }
};


class MCP: public Prox{
private:
    double gamma; // gamma_MCP >= 1
public:
    MCP(double g = 3){
        MoMALogger::debug("Initializing MCP proximal operator object");

        if(g < 1){
            MoMALogger::error("Non-convexity parameter (gamma) must be larger than or equal to 1 for MCP");
        }
        gamma = g;
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
        return z % sgn;
    }
};

template<class T>
class NonNegativeProx : public T{
public:
    NonNegativeProx<T>(): T() {
        MoMALogger::debug("Initializing non-negative prox");
    };

    NonNegativeProx<T>(double g): T(g) {
        MoMALogger::debug("Initializing non-negative prox");
    };

    arma::vec prox(const arma::vec &x, double l){
        return arma::max(T::prox(x, l), arma::zeros(x.n_elem));
    }
};

typedef NonNegativeProx<Lasso> NonNegativeLasso;
typedef NonNegativeProx<SCAD>  NonNegativeSCAD;
typedef NonNegativeProx<MCP>   NonNegativeMCP;

#endif
