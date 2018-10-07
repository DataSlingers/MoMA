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
        double i_lambda_u,      // regularization level
        double i_lambda_v,
        Rcpp::List i_prox_arg_list_u,
        Rcpp::List i_prox_arg_list_v,
        
        /*
         * smoothness - enforced through constraints
         */
        double i_alpha_u,               // Smoothing levels
        double i_alpha_v,
        const arma::mat &Omega_u,       // Smoothing matrices
        const arma::mat &Omega_v,
        
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

    // do parameter selection using nested BIC
    Rcpp::List select_nestedBIC( 
        const arma::vec &alpha_u,
        const arma::vec &alpha_v,
        const arma::vec &lambda_u,
        const arma::vec &lambda_v,
        int n_search);

    // deflate u * v.t() out of X by the amount of d
    int deflate(double d);

    // check convergence
    int check_cnvrg();

    // change penalty level
    int reset(double newlambda_u,double newlambda_v,
                double newalpha_u,double newalpha_v);
};

#endif
