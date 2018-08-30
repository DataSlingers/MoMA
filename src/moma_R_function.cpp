#include "moma.h"

// This function finds rank-k svd

// [[Rcpp::export]]
Rcpp::List cpp_sfpca(
    const arma::mat &X,    // We should not change any variable in R, so const ref
    const arma::mat &w_v,
    const arma::mat &w_u,
    const arma::mat &Omega_u, // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    std::string P_u,
    std::string P_v,
    double gamma_u,
    double gamma_v,
    double lambda2_u,
    double lambda2_v,
    bool ADMM_u,
    bool ADMM_v,
    bool acc_u,
    bool acc_v,
    double prox_eps_u,
    double prox_eps_v,
    bool nonneg_u,
    bool nonneg_v,
    arma::vec group_u,
    arma::vec group_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int k = 1){

    // WARNING: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
              /* sparsity */
              P_u,
              P_v,
              lambda_u(0),
              lambda_v(0),
              gamma_u,
              gamma_v,
              /* non-negativity */
              nonneg_u,
              nonneg_v,
              /* grouping */
              group_u,
              group_v,
              /* smoothness */
              Omega_u,
              Omega_v,
              alpha_u(0),
              alpha_v(0),
              /* sparse fused lasso*/
              lambda2_u,
              lambda2_v,
              /* unordered fusion*/
              w_u,
              w_v,
              ADMM_u,
              ADMM_v,
              acc_u,
              acc_v,
              prox_eps_u,
              prox_eps_v,
              /* algorithm parameters */
              EPS,
              MAX_ITER,
              EPS_inner,
              MAX_ITER_inner,
              solver);

    int n_lu = lambda_u.n_elem;
    int n_lv = lambda_v.n_elem;
    int n_au = alpha_u.n_elem;
    int n_av = alpha_v.n_elem;

    int n_more_than_one = int(n_lv > 1) + int(n_lu > 1) + int(n_au > 1) + int(n_av > 1);
    if(n_more_than_one > 0){
        MoMALogger::error("We don't allow a range of parameters in finding a rank-k svd.");
    }
    // store results
    arma::mat U(X.n_rows,k);
    arma::mat V(X.n_cols,k);
    arma::vec d(k);

    // find k PCs
    for(int i = 0; i < k; i++){
        problem.solve();
        U.col(i) = problem.u;
        V.col(i) = problem.v;
        d(i) = arma::as_scalar(problem.u.t() * problem.X * problem.v);
        // deflate X
        if(i < k-1){
            problem.deflate(d(i));
        }
    }
    return Rcpp::List::create(
                    Rcpp::Named("lambda_u") = lambda_u,
                    Rcpp::Named("lambda_v") = lambda_v,
                    Rcpp::Named("alpha_u") = alpha_u,
                    Rcpp::Named("alpha_v") = alpha_v,
                    Rcpp::Named("u") = U,
                    Rcpp::Named("v") = V,
                    Rcpp::Named("d") = d);
}

// This function solves a squence of lambda's and alpha's
// [[Rcpp::export]]
Rcpp::List cpp_sfpca_grid(
    const arma::mat &X,    // We should not change any variable in R, so const ref
    const arma::mat &w_v,
    const arma::mat &w_u,
    const arma::mat &Omega_u, // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    std::string P_u,
    std::string P_v,
    double gamma_u,
    double gamma_v,
    double lambda2_u,
    double lambda2_v,
    bool ADMM_u,
    bool ADMM_v,
    bool acc_u,
    bool acc_v,
    double prox_eps_u,
    double prox_eps_v,
    bool nonneg_u,
    bool nonneg_v,
    arma::vec group_u,
    arma::vec group_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int k = 1){

    // We only allow changing two parameters
    int n_lu = lambda_u.n_elem;
    int n_lv = lambda_v.n_elem;
    int n_au = alpha_u.n_elem;
    int n_av = alpha_v.n_elem;

    int n_more_than_one = int(n_lv > 1) + int(n_lu > 1) + int(n_au > 1) + int(n_av > 1);
    if(n_more_than_one > 2){
        MoMALogger::error("We only allow changing two parameters.");
    }

    if(n_lv == 0 || n_lu == 0 || n_au == 0 || n_av == 0){
        MoMALogger::error("Please specify all four parameters.");
    }

    int n_total = n_lv * n_lu * n_au * n_av;

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
              /* sparsity */
              P_u,
              P_v,
              lambda_u(0),
              lambda_v(0),
              gamma_u,
              gamma_v,
              /* non-negativity */
              nonneg_u,
              nonneg_v,
              /* grouping */
              group_u,
              group_v,
              /* smoothness */
              Omega_u,
              Omega_v,
              alpha_u(0),
              alpha_v(0),
              /* sparse fused lasso*/
              lambda2_u,
              lambda2_v,
              /* unordered fusion*/
              w_u,
              w_v,
              ADMM_u,
              ADMM_v,
              acc_u,
              acc_v,
              prox_eps_u,
              prox_eps_v,
              /* algorithm parameters */
              EPS,
              MAX_ITER,
              EPS_inner,
              MAX_ITER_inner,
              solver);

    // store results
    arma::mat U(X.n_rows,n_total);
    arma::mat V(X.n_cols,n_total);
    arma::vec d(n_total);

    int problem_id = 0;
    for(int i = 0; i < n_lu; i++){
        for(int j = 0; j < n_lv; j++){
            for(int k = 0; k < n_au; k++){
                for(int m = 0; m < n_av; m++){
                    MoMALogger::info("Setting up model:")
                                << " lambda_u " << lambda_u(i)
                                << " lambda_v " << lambda_v(j)
                                << " alpha_u " << alpha_u(k)
                                << " alpha_v " << alpha_v(m);

                    problem.reset(lambda_u(i),lambda_v(j),alpha_u(k),alpha_v(m));

                    // `solve` method use the result from last
                    // iteration as starting point
                    problem.solve();
                    U.col(problem_id) = problem.u;
                    V.col(problem_id) = problem.v;
                    d(problem_id) = arma::as_scalar(problem.u.t() * problem.X * problem.v);

                    problem_id++;
                }
            }
        }
    }
    if(problem_id != n_total){
        MoMALogger::error("Internal error: solution not found for all grid points.");
    }
    return Rcpp::List::create(
                    Rcpp::Named("lambda_u") = lambda_u,
                    Rcpp::Named("lambda_v") = lambda_v,
                    Rcpp::Named("alpha_u") = alpha_u,
                    Rcpp::Named("alpha_v") = alpha_v,
                    Rcpp::Named("u") = U,
                    Rcpp::Named("v") = V,
                    Rcpp::Named("d") = d);
}
