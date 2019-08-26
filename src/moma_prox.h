// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

#ifndef MOMA_PROX
#define MOMA_PROX 1

#include "moma_base.h"
#include "moma_logging.h"
#include "moma_prox_flsadp.h"
#include "moma_prox_fusion_util.h"
#include "moma_prox_sortedL1.h"

#define MAX(a, b) (a) > (b) ? (a) : (b)
#define MIN(a, b) (a) < (b) ? (a) : (b)
#define THRES_P(x, l) (MAX(x - l, 0.0))  // shrink a positive value by `l`

inline arma::vec soft_thres(const arma::vec &x, double l)
{
    return arma::sign(x) % arma::max(abs(x) - l, arma::zeros(arma::size(x)));
}

// soft-thresholding a non-negative vector
// all negative values will be set 0
inline arma::vec soft_thres_p(const arma::vec &x, double l)
{
    return arma::max(x - l, arma::zeros<arma::vec>(x.n_elem));
}

// An abstract class, member functions are implemeted in derived classes
class Prox
{
  public:
    virtual arma::vec operator()(const arma::vec &x, double l) = 0;
    virtual ~Prox()                                            = default;
    virtual int df(const arma::vec &x)                         = 0;
};

class NullProx : public Prox
{
  public:
    NullProx();
    arma::vec operator()(const arma::vec &x, double l);
    ~NullProx();
    int df(const arma::vec &x);
};

class Lasso : public Prox
{
  public:
    Lasso();
    arma::vec operator()(const arma::vec &x, double l);
    ~Lasso();
    int df(const arma::vec &x);
};

class SLOPE : public Prox
{
    arma::vec lambda;

  public:
    SLOPE(int dim);
    arma::vec operator()(const arma::vec &x, double l);
    ~SLOPE();
    int df(const arma::vec &x);
};

class NonNegativeLasso : public Prox
{
  public:
    NonNegativeLasso();
    arma::vec operator()(const arma::vec &x, double l);
    ~NonNegativeLasso();
    int df(const arma::vec &x);
};

class SCAD : public Prox
{
  protected:
    double gamma;  // gamma_SCAD >= 2
  public:
    SCAD(double g = 3.7);
    ~SCAD();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class NonNegativeSCAD : public SCAD
{
  public:
    NonNegativeSCAD(double g = 3.7);
    ~NonNegativeSCAD();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class MCP : public Prox
{
  protected:
    double gamma;  // gamma_MCP >= 1
  public:
    MCP(double g = 3);
    ~MCP();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class NonNegativeMCP : public MCP
{
  public:
    NonNegativeMCP(double g = 3);
    ~NonNegativeMCP();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class GrpLasso : public Prox
{
  protected:
    arma::vec group;
    int n_grp;  // number of gourps
  public:
    GrpLasso(const arma::vec &grp);
    ~GrpLasso();
    arma::vec operator()(const arma::vec &x, double l);
    arma::vec vec_prox(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class NonNegativeGrpLasso : public GrpLasso
{
  public:
    NonNegativeGrpLasso(const arma::vec &grp);
    ~NonNegativeGrpLasso();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class OrderedFusedLasso : public Prox
{
  public:
    OrderedFusedLasso();
    ~OrderedFusedLasso();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class OrderedFusedLassoDP : public OrderedFusedLasso
{
  public:
    OrderedFusedLassoDP();
    ~OrderedFusedLassoDP();
    arma::vec operator()(const arma::vec &x, double l);
};

class SparseFusedLasso : public Prox
{
  private:
    OrderedFusedLasso fg;
    double lambda2;  // lambda2 is the level of penalty on
                     // the absolute values of the coefficients

  public:
    SparseFusedLasso(double);
    ~SparseFusedLasso();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

class Fusion : public Prox
{
  private:
    arma::vec weight;
    bool ADMM;
    bool acc;
    double prox_eps;
    arma::vec start_point;

  public:
    Fusion(const arma::mat &input_w = arma::zeros<arma::mat>(0, 0),
           bool input_ADMM          = 1,
           bool input_acc           = 1,
           double input_prox_eps    = 1e-10);
    ~Fusion();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

// Its implementation is in `moma_prox_l1tf.cpp`
class L1TrendFiltering : public Prox
{
  private:
    int k;  // k \in 0,1,2, corresponding to fused lasso, linear tf and
            // third diff amt

    // The backtracking parameters
    // shrink stepsize by `bata`
    // if f(x + stepsize * dx) >= (1 - alpha * step) * f(x)
    // Ref:
    // http://www.stat.cmu.edu/~ryantibs/convexopt-F15/lectures/16-primal-dual.pdf
    // page 12
    static constexpr double alpha    = 0.01;
    static constexpr double beta     = 0.5;
    static const int MAX_ITER        = 200;
    static const int MAX_BT_ITER     = 5;
    static constexpr double prox_eps = 1e-10;

    arma::mat D;

  public:
    // n is the dim of the problem, k the degree of differences
    L1TrendFiltering(int n = -1, int i_k = 1);
    ~L1TrendFiltering();
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

// A handle class that deals with matching proximal operators
// and constructing and releasing the pointer
class ProxOp
{
  private:
    Prox *p;

  public:
    ProxOp() { p = nullptr; }

    ProxOp(Rcpp::List prox_arg_list, int dim);

    ~ProxOp() { delete p; }
    arma::vec operator()(const arma::vec &x, double l);
    int df(const arma::vec &x);
};

#endif
