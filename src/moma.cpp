// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma.h"

MoMA::MoMA(const arma::mat &i_X, // Pass X_ as a reference to avoid copy
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
    double i_alpha_u,           // Smoothing levels
    double i_alpha_v,
    const arma::mat &Omega_u,   // Smoothing matrices
    const arma::mat &Omega_v,

    /*
    * Algorithm parameters:
    */
    double i_EPS,
    long i_MAX_ITER,
    double i_EPS_inner,
    long i_MAX_ITER_inner,
    std::string i_solver):
    n(i_X.n_rows),
    p(i_X.n_cols),
    alpha_u(i_alpha_u),
    alpha_v(i_alpha_v),
    lambda_u(i_lambda_u),
    lambda_v(i_lambda_v),
    X(i_X),                                         // make our copy of the data
    MAX_ITER(i_MAX_ITER),
    EPS(i_EPS),
    solver_u(
            i_solver,
            alpha_u,Omega_u,
            lambda_u,i_prox_arg_list_u,
            i_EPS_inner,i_MAX_ITER_inner,i_X.n_rows),
    solver_v(
            i_solver,
            alpha_v,Omega_v,
            lambda_v,i_prox_arg_list_v,
            i_EPS_inner,i_MAX_ITER_inner,i_X.n_cols)
     // const reference must be passed to initializer list
{
    MoMALogger::info("Initializing MoMA object:")
    << " lambda_u " << lambda_u
    << " lambda_v " << lambda_v
    << " alpha_u " << alpha_u
    << " alpha_v " << alpha_v
    << " P_u " << Rcpp::as<std::string>(i_prox_arg_list_u["P"])
    << " P_v " << Rcpp::as<std::string>(i_prox_arg_list_v["P"]);
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
        MoMALogger::error("Cannot deflate by non-positive factor.");
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
        MoMALogger::debug("Outer loop No.") << iter << "--" << tol;
    }
    
    MoMALogger::info("--Finish iter: ") << iter << "---" ;
    check_cnvrg();
}

int MoMA::check_cnvrg(){
    if(iter >= MAX_ITER){
        MoMALogger::warning("No convergence in MoMA!") 
            << " lambda_u " << lambda_u
            << " lambda_v " << lambda_v
            << " alpha_u " << alpha_u
            << " alpha_v " << alpha_v;
    }
    return 0;
} 

int MoMA::reset(double newlambda_u,double newlambda_v,
                double newalpha_u,double newalpha_v){

    solver_u.reset(newlambda_u,newalpha_u);
    solver_v.reset(newlambda_v,newalpha_v);
    return 0;
}
