// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#ifndef MOMA_H
#define MOMA_H 1
// Global #includes and #defines
#include "moma_base.h"

// Logging code
#include "moma_logging.h"

// Prox operators
#include "moma_solver.h"

// Prototypes
// moma_logging.cpp
void moma_set_logger_level_cpp(int);
int moma_get_logger_level_cpp();
void moma_log_cpp(int, Rcpp::StringVector);

class MoMA{

private:
    /* matrix size */
    int n;      // rows
    int p;      // columns
    double alpha_u;
    double alpha_v;
    double lambda_u;
    double lambda_v;

public:
    // Our own copy of the data matrix
    // Modified when finding rank-k svd
    arma::mat X;
    // user-specified precisions
    int MAX_ITER;
    double EPS;
    // working precision and iterations
    int iter;
    double tol;
    // Results -- will be modified during iterations and copied back to R
    arma::vec u;
    arma::vec v;

    PR_solver solver_u;
    PR_solver solver_v;
 
    // Parse user input into a MoMA object which defines the problem and algorithm
    // used to solve it.
    //
    // TODO: Decouple problem defintion and algorithmic choices
    //
    MoMA(const arma::mat &X_,   // Pass X_ as a reference to avoid copy
        /*
         * sparsity - enforced through penalties
         */
        std::string P_u,        // Sparsity penalty info
        std::string P_v,
        double i_lambda_u,      // regularization level
        double i_lambda_v,
        double gamma_u,    
        double gamma_v,         // Non-convexity parameter
        bool nonneg_u,          // Non-negativity indicator
        bool nonneg_v,
        /*
        * grouping
        */
        const arma::vec &group_u,
        const arma::vec &group_v,
        /*
         * smoothness - enforced through constraints
         */
        const arma::mat &Omega_u,       // Smoothing matrices
        const arma::mat &Omega_v,
        double i_alpha_u,               // Smoothing levels
        double i_alpha_v,
        /*
         * unordered fusion
         */
        const arma::mat &w_u,
        const arma::mat &w_v,
        bool ADMM_u,
        bool ADMM_v,
        bool acc_u,
        bool acc_v,
        double prox_eps_u,
        double prox_eps_v,
        /*
         * Algorithm parameters:
         */
        double i_EPS,
        long i_MAX_ITER,
        double i_EPS_inner,
        long i_MAX_ITER_inner,
        std::string i_solver);

    // solve sfpca by iteratively solving
    // penalized regressions
    void solve();

    // defalte u * v.t() out of X by the amount of d
    int deflate(double d);

    // check convergence
    int check_cnvrg();
};

#endif
