// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-

#ifndef MOMA_SOLVER
#define MOMA_SOLVER 1

#include "moma_base.h"
#include "moma_logging.h"
#include "moma_prox.h"

// Penalized regression solver
// min_u || y - u || + lambda * P(u) s.t. || u ||_S <= 1
// S = I + alpha * Omega
class _PR_solver
{
  protected:
    int dim;  // dimension of the PR problem
    double lambda;
    double alpha;
    double L;
    const arma::mat &Omega;
    // S = I + alpha * Omega for u, v smoothing
    arma::mat S;
    bool is_S_idmat;  // indicator of alpha == 0.0 <=> S == I

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
    arma::vec g(const arma::vec &v,
                const arma::vec &y,
                double step_size,
                const arma::mat &S,
                bool is_S_idmat);
    arma::vec normalize(const arma::vec &u);

    // user-specified precision and max iterations
    double EPS;
    int MAX_ITER;

  public:
    explicit _PR_solver(
        // smoothness
        double i_alpha,
        const arma::mat &i_Omega,
        // sparsity
        double i_lambda,
        Rcpp::List prox_arg_list,
        // algorithm settings
        double i_EPS,
        int i_MAX_ITER,
        int i_dim);

    // Used when solving for a bunch of lambda's and alpha's
    int set_penalty(double new_lambda, double new_alpha);
    double bic(arma::vec y, const arma::vec &est);
    virtual ~_PR_solver()                                              = default;
    virtual arma::vec solve(arma::vec y, const arma::vec &start_point) = 0;
    void check_convergence(int iter, double tol);
};

class ISTA : public _PR_solver
{
  public:
    ISTA(double i_alpha,
         const arma::mat &i_Omega,
         double i_lambda,
         Rcpp::List prox_arg_list,
         double i_EPS,
         int i_MAX_ITER,
         int dim)
        : _PR_solver(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim)
    {
        MoMALogger::debug("Initializing a ISTA solver.");
    };
    arma::vec solve(arma::vec y, const arma::vec &start_point);
    ~ISTA() { MoMALogger::debug("Releasing a ISTA object"); }
};

class FISTA : public _PR_solver
{
  public:
    FISTA(double i_alpha,
          const arma::mat &i_Omega,
          double i_lambda,
          Rcpp::List prox_arg_list,
          double i_EPS,
          int i_MAX_ITER,
          int dim)
        : _PR_solver(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim)
    {
        MoMALogger::debug("Initializing a FISTA solver.");
    };
    arma::vec solve(arma::vec y, const arma::vec &start_point);
    ~FISTA() { MoMALogger::debug("Releasing a FISTA object"); }
};

class OneStepISTA : public _PR_solver
{
  public:
    OneStepISTA(double i_alpha,
                const arma::mat &i_Omega,
                double i_lambda,
                Rcpp::List prox_arg_list,
                double i_EPS,
                int i_MAX_ITER,
                int dim)
        : _PR_solver(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim)
    {
        MoMALogger::debug("Initializing an one-step ISTA solver.");
    };
    arma::vec solve(arma::vec y, const arma::vec &start_point);
    ~OneStepISTA() { MoMALogger::debug("Releasing a OneStepISTA object"); }
};

// A handle class
class PR_solver
{
  private:
    _PR_solver *prs;

  public:
    PR_solver(
        // a string saying which algorithm to use
        const std::string &algorithm_string,
        // same as class _PR_solver
        double i_alpha,
        const arma::mat &i_Omega,
        double i_lambda,
        Rcpp::List prox_arg_list,
        double i_EPS,
        int i_MAX_ITER,
        int dim);

    // wrap operations in _PR_solver class
    arma::vec solve(arma::vec y, const arma::vec &start_point);
    double bic(arma::vec y, const arma::vec &est);
    int set_penalty(double new_lambda, double new_alpha);

    ~PR_solver() { delete prs; }
};

#endif
