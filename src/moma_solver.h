// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef MOMA_SOLVER
#define MOMA_SOLVER 1

#include "moma_base.h"
#include "moma_logging.h"
#include "moma_prox.h"

// Penalized regression solver
// min_u || y - u || + lambda * P(u) s.t. || u ||_S <= 1
// S = I + alpha * Omega
class _PR_solver{

protected:
    int dim;                    // dimension of the PR problem
    double lambda;
    double alpha;
    double L;
    const arma::mat &Omega;
     // S = I + alpha * Omega for u, v smoothing
    arma::mat S;
    bool I; // indicator of alpha == 0.0 <=> S == I


    // Step size for proximal gradient algorithm
    //   - since this is a linear model internally, we can used a fixed
    //     step size without backtracking
    double grad_step_size;
    double prox_step_size;
    // A proximal operator for sparsity inducing penalties
    //
    // Note that currently the threshold level is not defined in the Prox object
    ProxOp p;
    // A gradient operator
    arma::vec g(
        const arma::vec &v, const arma::vec &y,
        double step_size, const arma::mat &S, bool I);
    arma::vec normalize(const arma::vec &u);


    // user-specified precision and max iterations
    double EPS;
    int MAX_ITER;
    // working precision and max iterations
    double tol;
    int iter;
    int check_cnvrg();

public:
    explicit _PR_solver(
        // smoothness
        double i_alpha, const arma::mat &i_Omega,
        // sparsity
        double i_lambda, const std::string &sparsity_string, double gamma,
        const arma::vec &group, const arma::mat &w, bool ADMM,
        bool acc, double prox_eps, bool nonneg,
        // algorithm settings
        double i_EPS, int i_MAX_ITER, int i_dim);

    // Used when solving for a bunch of lambda's and alpha's
    int reset(double new_lambda, double new_alpha);
    virtual ~_PR_solver() = default;
    virtual arma::vec solve(arma::vec y, const arma::vec &start_point) = 0;
};

class ISTA: public _PR_solver{
public:
    ISTA(
        double i_alpha, const arma::mat &i_Omega, double i_lambda,
        const std::string &sparsity_string, double gamma, const arma::vec &group,
        const arma::mat &w, bool ADMM, bool acc, double prox_eps, bool nonneg,
        double i_EPS, int i_MAX_ITER, int dim)
        : _PR_solver(
                i_alpha,i_Omega,i_lambda,sparsity_string,gamma,
                group,w,ADMM,acc,prox_eps,nonneg,i_EPS,i_MAX_ITER,dim)
    {
        MoMALogger::debug("Initializing a ISTA solver.");
    };
    arma::vec solve(arma::vec y,const arma::vec &start_point);
};

class FISTA: public _PR_solver{
public:
    FISTA(
        double i_alpha, const arma::mat &i_Omega, double i_lambda,
        const std::string &sparsity_string, double gamma, const arma::vec &group,
        const arma::mat &w, bool ADMM, bool acc, double prox_eps, bool nonneg,
        double i_EPS, int i_MAX_ITER, int dim)
        :_PR_solver(
                i_alpha,i_Omega,i_lambda,sparsity_string,gamma,
                group,w,ADMM,acc,prox_eps,nonneg,i_EPS,i_MAX_ITER,dim)
    {
        MoMALogger::debug("Initializing a FISTA solver.");
    };
    arma::vec solve(arma::vec y,const arma::vec &start_point);
};

class OneStepISTA: public _PR_solver{
public:
    OneStepISTA(
                double i_alpha, const arma::mat &i_Omega, double i_lambda,
                const std::string &sparsity_string, double gamma, const arma::vec &group,
                const arma::mat &w, bool ADMM, bool acc, double prox_eps, bool nonneg,
                double i_EPS, int i_MAX_ITER, int dim)
        :_PR_solver(
                i_alpha,i_Omega,i_lambda,sparsity_string,gamma,
                group,w,ADMM,acc,prox_eps,nonneg,i_EPS,i_MAX_ITER,dim)
    {
        MoMALogger::debug("Initializing an one step ISTA solver.");
    };
    arma::vec solve(arma::vec y,const arma::vec &start_point);
};

// A handle class
class PR_solver{
private:
    _PR_solver *prs;
public:
    PR_solver(
        // a string saying which algorithm to use
        const std::string &algorithm_string,
        // same as class _PR_solver
        double i_alpha, const arma::mat &i_Omega,
        double i_lambda, const std::string &sparsity_string, double gamma,
        const arma::vec &group, const arma::mat &w, bool ADMM, 
        bool acc, double prox_eps, bool nonneg,
        double i_EPS, int i_MAX_ITER, int dim);

    // wrap operations in _PR_solver class
    arma::vec solve(arma::vec y, const arma::vec &start_point);
    int reset(double new_lambda, double new_alpha);

    ~PR_solver(){
        delete prs;
    }
};

#endif
