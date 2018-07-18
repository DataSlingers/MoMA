// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef MOMA_PROX
#define MOMA_PROX 1

#include "moma_base.h"
#include "moma_logging.h"
#include "moma_prox_fusion_util.h"
#define MAX(a,b) (a)>(b)?(a):(b)
#define THRES_P(x,l) (MAX(x-l,0.0)) // shrink a positive value by `l`

inline arma::vec soft_thres(const arma::vec &x, double l){
    return arma::sign(x) % arma::max(abs(x) - l, arma::zeros(arma::size(x)));
}

// soft-thresholding a non-negative vector
// all negative values will be set 0
inline arma::vec soft_thres_p(const arma::vec &x, double l){
    return arma::max(x - l, arma::zeros<arma::vec>(x.n_elem));
}

// An abstract class, member functions are implemeted in derived classes
class Prox{
public:
    virtual arma::vec operator()(const arma::vec &x, double l) = 0;
    virtual ~Prox() = default;
};

class NullProx: public Prox{
public:
    NullProx();
    virtual arma::vec operator()(const arma::vec &x, double l);
    virtual ~NullProx();
};

class Lasso: public Prox{
public:
    Lasso();
    arma::vec operator()(const arma::vec &x, double l);
    ~Lasso();
};

class NonNegativeLasso: public Prox{
public:
    NonNegativeLasso();
    arma::vec operator()(const arma::vec &x, double l);
    ~NonNegativeLasso();
};

class SCAD: public Prox{
protected:
    double gamma; // gamma_SCAD >= 2
public:
    SCAD(double g = 3.7);
    ~SCAD();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
};

class NonNegativeSCAD: public SCAD{
public:
    NonNegativeSCAD(double g = 3.7);
    ~NonNegativeSCAD();
    arma::vec operator()(const arma::vec &x, double l);
};

class MCP: public Prox{
protected:
    double gamma; // gamma_MCP >= 1
public:
    MCP(double g = 3);
    ~MCP();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
};

class NonNegativeMCP: public MCP{
public:
    NonNegativeMCP(double g = 3);
    ~NonNegativeMCP();
    arma::vec operator()(const arma::vec &x, double l);
};

class GrpLasso: public Prox{
protected:
    arma::vec group;
    int n_grp; // number of gourps
    arma::umat D;  // Probably not using sparse matrix would be faster, TODO
                    // a boolean matrix, D \in R^{g \times p}, g is the number of groups, p the number of parameters.
                    // D_ji = 1 means \beta_i in group j.
                    // should be integer, probably use arma::sp_umat; it will cause error though, when it multipies a vec
public:
    GrpLasso(const arma::vec &grp);
    ~GrpLasso();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
};

class NonNegativeGrpLasso: public GrpLasso{
public:
    NonNegativeGrpLasso(const arma::vec &grp);
    ~NonNegativeGrpLasso();
    arma::vec operator()(const arma::vec &x, double l);
};

class OrderedFusion: public Prox{
public:
    OrderedFusion();
    ~OrderedFusion();
    arma::vec operator()(const arma::vec &x, double l);       
};

#endif
