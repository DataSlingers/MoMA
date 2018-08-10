#include "moma.h"

// This function finds rank-k svd

// [[Rcpp::export]]
Rcpp::List cpp_sfpca(
    const arma::mat &X,    // We should not change any variable in R, so const ref
    const arma::mat &w_v,
    const arma::mat &w_u,
    const arma::mat &Omega_u, // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    double alpha_u,
    double alpha_v,
    double lambda_u,
    double lambda_v,
    std::string P_u,
    std::string P_v,
    double gamma,
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
    std::string solver,
    int k = 1){

    MoMA problem(X,
              /* sparsity */
              P_v,
              P_u,
              lambda_v,
              lambda_u,
              gamma,
              /* non-negativity */
              nonneg_u,
              nonneg_v,
              /* grouping */
              group_u,
              group_v,
              /* smoothness */
              Omega_u,
              Omega_v,
              alpha_u,
              alpha_v,
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
              solver);

    // store results
    arma::mat U(X.n_rows,k);
    arma::mat V(X.n_cols,k);
    arma::vec d(k);

    // find k PCs
    for(int i = 0; i < k; i++){
        problem.solve();
        U.col(i) = problem.u;
        V.col(i) = problem.v;
        d(i) = arma::as_scalar(trans(problem.u) * problem.X * problem.v);
        // deflate X
        if(i < k-1){
            problem.deflate(d(i));
        }
    }
    return Rcpp::List::create(
                    Rcpp::Named("u") = U,
                    Rcpp::Named("v") = V,
                    Rcpp::Named("d") = d);
}

// This function solves a squence of lambda
// TODO