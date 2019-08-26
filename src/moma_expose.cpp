#include "moma.h"

// This files expose four member functions of MoMA
// 1. MoMA::multi_rank (see function `cpp_moma_multi_rank`)
// 2. MoMA::grid_search (see function `cpp_moma_grid_search`)
// 3. MoMA::criterion_search (see function `cpp_moma_criterion_search`)

// [[Rcpp::export]]
Rcpp::List cpp_moma_multi_rank(
    const arma::mat &X,  // We should not change any variable in R, so const ref
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::mat &Omega_u,  // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    const Rcpp::List &prox_arg_list_u,
    const Rcpp::List &prox_arg_list_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int rank = 1)
{
    // WARNING: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
                 /* sparsity */
                 lambda_u(0), lambda_v(0), prox_arg_list_u, prox_arg_list_v,
                 /* smoothness */
                 alpha_u(0), alpha_v(0), Omega_u, Omega_v,
                 /* algorithm parameters */
                 EPS, MAX_ITER, EPS_inner, MAX_ITER_inner, solver);

    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;

    int n_more_than_one =
        int(n_lambda_v > 1) + int(n_lambda_u > 1) + int(n_alpha_u > 1) + int(n_alpha_v > 1);
    if (n_more_than_one > 0)
    {
        MoMALogger::error("We don't allow a range of parameters in finding a rank-k svd.");
    }
    return problem.multi_rank(rank, problem.u, problem.v);
}

// This function solves a squence of lambda's and alpha's
// [[Rcpp::export]]
Rcpp::List cpp_moma_grid_search(
    const arma::mat &X,  // We should not change any variable in R, so const ref
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::mat &Omega_u,  // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    const Rcpp::List &prox_arg_list_u,
    const Rcpp::List &prox_arg_list_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int rank = 1)  // `rank` is not used
{
    if (rank != 1)
    {
        MoMALogger::error("Then `rank` argument in `cpp_moma_multi_rank` should not be specified.");
    }
    // We only allow changing two parameters
    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;

    int n_more_than_one =
        int(n_lambda_v > 1) + int(n_lambda_u > 1) + int(n_alpha_u > 1) + int(n_alpha_v > 1);
    if (n_more_than_one > 2)
    {
        MoMALogger::error("We only allow changing two parameters.");
    }

    if (n_lambda_v == 0 || n_lambda_u == 0 || n_alpha_u == 0 || n_alpha_v == 0)
    {
        MoMALogger::error("Please specify all four parameters.");
    }

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
                 /* sparsity */
                 lambda_u(0), lambda_v(0), prox_arg_list_u, prox_arg_list_v,
                 /* smoothness */
                 alpha_u(0), alpha_v(0), Omega_u, Omega_v,
                 /* algorithm parameters */
                 EPS, MAX_ITER, EPS_inner, MAX_ITER_inner, solver);

    // store results
    return problem.grid_search(alpha_u, lambda_u, alpha_v, lambda_v, problem.u, problem.v);
}

// This function solves a squence of lambda's and alpha's
// [[Rcpp::export]]
Rcpp::List cpp_moma_criterion_search(
    const arma::mat &X,  // We should not change any variable in R, so const ref
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::mat &Omega_u,  // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    const Rcpp::List &prox_arg_list_u,
    const Rcpp::List &prox_arg_list_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int rank = 1)  // rank not used
{
    if (rank != 1)
    {
        MoMALogger::error("Then `rank` argument in `cpp_moma_multi_rank` should not be specified.");
    }
    // We only allow changing two parameters
    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;

    int n_more_than_one =
        int(n_lambda_v > 1) + int(n_lambda_u > 1) + int(n_alpha_u > 1) + int(n_alpha_v > 1);
    if (n_more_than_one > 2)
    {
        MoMALogger::error("We only allow changing two parameters.");
    }

    if (n_lambda_v == 0 || n_lambda_u == 0 || n_alpha_u == 0 || n_alpha_v == 0)
    {
        MoMALogger::error("Please specify all four parameters.");
    }

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
                 /* sparsity */
                 lambda_u(0), lambda_v(0), prox_arg_list_u, prox_arg_list_v,
                 /* smoothness */
                 alpha_u(0), alpha_v(0), Omega_u, Omega_v,
                 /* algorithm parameters */
                 EPS, MAX_ITER, EPS_inner, MAX_ITER_inner, solver);

    return problem.criterion_search(alpha_u, lambda_u, alpha_v, lambda_v, problem.u, problem.v,
                                    EPS);
}

// This function solves a squence of lambda's and alpha's
// [[Rcpp::export]]
Rcpp::List cpp_multirank_BIC_grid_search(
    const arma::mat &X,  // We should not change any variable in R, so const ref
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::mat &Omega_u,  // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    const Rcpp::List &prox_arg_list_u,
    const Rcpp::List &prox_arg_list_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int deflation_scheme       = 1,
    int select_scheme_alpha_u  = 0,  // 0 = grid, 1 = bic
    int select_scheme_alpha_v  = 0,
    int select_scheme_lambda_u = 0,
    int select_scheme_lambda_v = 0,
    int max_bic_iter           = 5,
    int rank                   = 1)
{
    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;

    if (n_lambda_v == 0 || n_lambda_u == 0 || n_alpha_u == 0 || n_alpha_v == 0)
    {
        MoMALogger::error("Please specify all four parameters.");
    }

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
                 /* sparsity */
                 lambda_u(0), lambda_v(0), prox_arg_list_u, prox_arg_list_v,
                 /* smoothness */
                 alpha_u(0), alpha_v(0), Omega_u, Omega_v,
                 /* algorithm parameters */
                 EPS, MAX_ITER, EPS_inner, MAX_ITER_inner, solver,
                 static_cast<DeflationScheme>(deflation_scheme));

    return problem.grid_BIC_mix(
        alpha_u, alpha_v, lambda_u, lambda_v, static_cast<SelectionScheme>(select_scheme_alpha_u),
        static_cast<SelectionScheme>(select_scheme_alpha_v),
        static_cast<SelectionScheme>(select_scheme_lambda_u),
        static_cast<SelectionScheme>(select_scheme_lambda_v), max_bic_iter, rank);
}

// [[Rcpp::export]]
Rcpp::List cca(const arma::mat &X,  // We should not change any variable in R, so const ref
               const arma::mat &Y,
               const arma::vec &alpha_u,
               const arma::vec &alpha_v,
               const arma::mat &Omega_u,  // Default values for these matrices should be set in R
               const arma::mat &Omega_v,
               const arma::vec &lambda_u,
               const arma::vec &lambda_v,
               const Rcpp::List &prox_arg_list_u,
               const Rcpp::List &prox_arg_list_v,
               double EPS,
               long MAX_ITER,
               double EPS_inner,
               long MAX_ITER_inner,
               std::string solver,
               int deflation_scheme,            // refer to `DeflationScheme` in moma_base.h
               int select_scheme_alpha_u  = 0,  // 0 means grid, 1 means BIC search
               int select_scheme_alpha_v  = 0,
               int select_scheme_lambda_u = 0,
               int select_scheme_lambda_v = 0,
               int max_bic_iter           = 5,
               int rank                   = 1)
{
    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;

    if (n_lambda_v == 0 || n_lambda_u == 0 || n_alpha_u == 0 || n_alpha_v == 0)
    {
        MoMALogger::error("Please specify all four parameters.");
    }

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X, Y,
                 /* sparsity */
                 lambda_u(0), lambda_v(0), prox_arg_list_u, prox_arg_list_v,
                 /* smoothness */
                 alpha_u(0), alpha_v(0), Omega_u, Omega_v,
                 /* algorithm parameters */
                 EPS, MAX_ITER, EPS_inner, MAX_ITER_inner, solver,
                 static_cast<DeflationScheme>(deflation_scheme));

    return problem.grid_BIC_mix(
        alpha_u, alpha_v, lambda_u, lambda_v, static_cast<SelectionScheme>(select_scheme_alpha_u),
        static_cast<SelectionScheme>(select_scheme_alpha_v),
        static_cast<SelectionScheme>(select_scheme_lambda_u),
        static_cast<SelectionScheme>(select_scheme_lambda_v), max_bic_iter, rank);
}
