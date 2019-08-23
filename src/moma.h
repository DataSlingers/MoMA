// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
#ifndef MOMA_H
#define MOMA_H 1
// Global #includes and #defines
#include "moma_base.h"

// Logging code
#include "moma_logging.h"

// Prox operators
#include "moma_solver.h"

// BIC searcher
#include "moma_solver_BICsearch.h"

// 4-D list
#include "moma_fivedlist.h"

// Prototypes
// moma_logging.cpp
void moma_set_logger_level_cpp(int);
int moma_get_logger_level_cpp();
void moma_log_cpp(int, Rcpp::StringVector);

class MoMA
{
  private:
    /* matrix size */
    int n;  // rows
    int p;  // columns
    double alpha_u;
    double alpha_v;
    double lambda_u;
    double lambda_v;
    bool is_initialzied;  // only MoMA::initialze_uv() sets it to true
    bool is_solved;       // only MoMA::solve() sets it true.
    arma::mat X;          // on X we perform the algorithm

    arma::mat X_working;  // keep track of deflated matrices
    arma::mat Y_working;

    arma::mat X_original;  // const
    arma::mat Y_original;

    const DeflationScheme ds;

    // they are used in MoMA::evaluate_loss
    arma::mat Omega_u;
    arma::mat Omega_v;

  public:
    // Receiver a grid of parameters
    // and perform greedy BIC search. Initial points
    // must be specified.
    BIC_searcher bicsr_u;
    BIC_searcher bicsr_v;

    // Our own copy of the data matrix
    // Modified when finding rank-k svd
    // user-specified precisions
    int MAX_ITER;
    double EPS;
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
    MoMA(const arma::mat &X_,  // Pass X_ as a reference to avoid copy
         /*
          * sparsity - enforced through penalties
          */
         double i_lambda_u,  // regularization level
         double i_lambda_v,
         Rcpp::List i_prox_arg_list_u,
         Rcpp::List i_prox_arg_list_v,

         /*
          * smoothness - enforced through constraints
          */
         double i_alpha_u,  // Smoothing levels
         double i_alpha_v,
         const arma::mat &i_Omega_u,  // Smoothing matrices
         const arma::mat &i_Omega_v,

         /*
          * Algorithm parameters:
          */
         double i_EPS,
         long i_MAX_ITER,
         double i_EPS_inner,
         long i_MAX_ITER_inner,
         std::string i_solver,
         DeflationScheme i_ds = DeflationScheme::PCA_Hotelling);

    MoMA(
        // Pass X_ as a reference to avoid copy
        const arma::mat &X_working,
        const arma::mat &Y_working,
        /*
         * sparsity - enforced through penalties
         */
        double i_lambda_u,  // regularization level
        double i_lambda_v,
        Rcpp::List i_prox_arg_list_u,
        Rcpp::List i_prox_arg_list_v,

        /*
         * smoothness - enforced through constraints
         */
        double i_alpha_u,  // Smoothing levels
        double i_alpha_v,
        const arma::mat &i_Omega_u,  // Smoothing matrices
        const arma::mat &i_Omega_v,

        /*
         * Algorithm parameters:
         */
        double i_EPS,
        long i_MAX_ITER,
        double i_EPS_inner,
        long i_MAX_ITER_inner,
        std::string i_solver,
        DeflationScheme i_ds);

    // solve sfpca by iteratively solving
    // penalized regressions
    void solve();

    double evaluate_loss();

    // Deflation happens in place, so MoMA::X is contaminated
    int deflate();

    int initialize_uv();

    // check convergence
    int check_convergence(int iter, double tol);

    // change penalty level
    int set_penalty(double newlambda_u, double newlambda_v, double newalpha_u, double newalpha_v);
    int reset_X();

    // following functions are implemented in `moma_level1.cpp`
    Rcpp::List criterion_search(const arma::vec &bic_au_grid,
                                const arma::vec &bic_lu_grid,
                                const arma::vec &bic_av_grid,
                                const arma::vec &bic_lv_grid,
                                arma::vec initial_u,
                                arma::vec initial_v,
                                double EPS_bic   = 1e-7,  // not very useful
                                int max_bic_iter = 5,
                                bool final_run   = true);

    Rcpp::List multi_rank(int rank, arma::vec initial_u, arma::vec initial_v);

    Rcpp::List grid_search(const arma::vec &alpha_u,
                           const arma::vec &lambda_u,
                           const arma::vec &alpha_v,
                           const arma::vec &lambda_v,
                           arma::vec initial_u,
                           arma::vec initial_v);

    Rcpp::List grid_BIC_mix(const arma::vec &alpha_u,
                            const arma::vec &alpha_v,
                            const arma::vec &lambda_u,
                            const arma::vec &lambda_v,
                            SelectionScheme select_scheme_alpha_u,
                            SelectionScheme select_scheme_alpha_v,
                            SelectionScheme select_scheme_lambda_u,
                            SelectionScheme select_scheme_lambda_v,
                            int max_bic_iter = 5,
                            int rank         = 1);
};

#endif
