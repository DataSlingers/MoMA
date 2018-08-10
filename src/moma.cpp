// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma.h"

MoMA::MoMA(const arma::mat &i_X, // Pass X_ as a reference to avoid copy
    /*
    * sparsity - enforced through penalties
    */
    std::string P_v,    // Sparsity penalty info
    std::string P_u,
    double i_lambda_v,  // regularization level
    double i_lambda_u,
    double gamma,       // Non-convexity parameter
    bool nonneg_u,      // Non-negativity indicator
    bool nonneg_v,
    /*
    * grouping
    */
    const arma::vec &group_u,
    const arma::vec &group_v,
    /*
    * smoothness - enforced through constraints
    */
    const arma::mat &Omega_u,   // Smoothing matrices
    const arma::mat &Omega_v,
    double i_alpha_u,           // Smoothing levels
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
    arma::uword i_MAX_ITER,
    std::string i_solver):

    alpha_u(i_alpha_u),
    alpha_v(i_alpha_v),
    lambda_u(i_lambda_u),
    lambda_v(i_lambda_v),
    X(i_X),
    MAX_ITER(i_MAX_ITER),
    EPS(i_EPS),
    solver_u(i_solver,alpha_u,Omega_u,lambda_u,P_u,gamma,group_u,w_u,ADMM_u,acc_u,prox_eps_u,nonneg_u,i_EPS,i_MAX_ITER),
    solver_v(i_solver,alpha_v,Omega_v,lambda_v,P_v,gamma,group_v,w_v,ADMM_v,acc_v,prox_eps_v,nonneg_v,i_EPS,i_MAX_ITER)
     // const reference must be passed to initializer list
{
    MoMALogger::info("Setting up model");

    // Step 2: Initialize to leading singular vectors
    //
    //         MoMA is a regularized SVD, which is a non-convex (bi-convex)
    //         problem, so we need to be cautious about initialization to
    //         avoid local-minima. Initialization at the SVD (global solution
    //         to the non-regularized problem) seems to be a good trade-off:
    //         for problems with little regularization, the MoMA solution will
    //         lie near the SVD solution; for problems with significant regularization
    //         the problem becomes more well-behaved and less sensitive to
    //         initialization
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    v = V.col(0);
    u = U.col(0);
};

int MoMA::deflate(double d){
    MoMALogger::warning("Deflating.");
    if(d <= 0.0){
        MoMALogger::warning("Defalte by non-positive number.");
    }
    X = X - d * u * v.t();
    // Re-initialize u and v after deflation
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    v = V.col(0);
    u = U.col(0);
    return d;
}

void MoMA::solve(){
    tol = 1;
    iter = 0;
    arma::vec oldu;
    arma::vec oldv;
    while(tol > EPS && iter < MAX_ITER){
        iter++;
        oldu = u;
        oldv = v;

        u = solver_u.solve(X*v, u);
        v = solver_v.solve(X.t()*u, v);

        tol = norm(oldu - u) / norm(oldu) + norm(oldv - v) / norm(oldv);
        MoMALogger::debug("Outer loop No.") << iter << "--"<< "%change " << tol;
    }
    MoMALogger::info("--Finish iter: ") << iter << "---" ;
}
