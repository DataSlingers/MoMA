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

    bicsr_u.bind(&solver_u, &PR_solver::bic);
    bicsr_v.bind(&solver_v, &PR_solver::bic);

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

Rcpp::List MoMA::select_nestedBIC( 
        const arma::vec &alpha_u,
        const arma::vec &alpha_v,
        const arma::vec &lambda_u,
        const arma::vec &lambda_v,
        int max_bic_iter = 5){          // suggested in the sfpca_nested_bic.m
    
    MoMALogger::info("Running nested BIC parameter selection.");
    tol = 1;
    iter = 0;
    arma::vec oldu;
    arma::vec oldv;
    arma::vec working_u = u;
    arma::vec working_v = v;
    double working_bic_u;
    double working_bic_v;
    double minbic_u;
    double minbic_v;
    double opt_alpha_u;
    double opt_alpha_v;
    double opt_lambda_u;
    double opt_lambda_v;

    while(tol > EPS && iter < MAX_ITER && iter < max_bic_iter){
        iter++;
        oldu = u;
        oldv = v;
        minbic_u = 1e+10;
        minbic_v = 1e+10;

        // choose lambda/alpha_u
        for(int i = 0; i < alpha_u.n_elem; i++){
            for(int j = 0; j < lambda_u.n_elem; j++){
                // Put lambda_u in the inner loop to avoid reconstructing S many times
                solver_u.reset(lambda_u(j),alpha_u(i));
                working_u     = solver_u.solve(X*v, working_u);
                working_bic_u = solver_u.bic(X*v, working_u);
                MoMALogger::debug("(now,min,la,al) = (") << working_bic_u << "," << minbic_u << "," << lambda_u(j) << "," << alpha_u(i) << ")";
                if(working_bic_u < minbic_u){
                    minbic_u     = working_bic_u;
                    u            = working_u;
                    opt_lambda_u = lambda_u(j);
                    opt_alpha_u  = alpha_u(i);
                }
            }
        }
        MoMALogger::message("Search No.") << iter << ", BIC(u) = " << minbic_u << 
                                        ", (al,lam) = (" << opt_alpha_u << 
                                        ", " << opt_lambda_u << ").";
        // choose lambda/alpha_v
        for(int i = 0; i < alpha_v.n_elem; i++){
            for(int j = 0; j < lambda_v.n_elem; j++){
                // Put lambda_v in the inner loop to avoid reconstructing S many times
                solver_v.reset(lambda_v(j),alpha_v(i));
                working_v     = solver_v.solve(X.t()*u, working_v);
                working_bic_v = solver_v.bic(X.t()*u, working_v);
                MoMALogger::debug("(now,min) = (") << working_bic_v << "," << minbic_v << "," << lambda_v(j) << "," << alpha_v(i) << ")";
                if(working_bic_v < minbic_v){
                    minbic_v     = working_bic_v;
                    v            = working_v;
                    opt_lambda_v = lambda_v(j);
                    opt_alpha_v  = alpha_v(i);
                }
            }
        }
        MoMALogger::message("Search No.") << iter << ", BIC(v) = " << minbic_v << 
                                                ", (al,lam) = (" << opt_alpha_v << 
                                                ", " << opt_lambda_v << ").";

        tol = norm(oldu - u) / norm(oldu) + norm(oldv - v) / norm(oldv);
        MoMALogger::debug("Outer loop No.") << iter << "--" << tol;
    }
    
    // A final run on the chosen set of parameters
    reset(opt_lambda_u,opt_lambda_v,opt_alpha_u,opt_alpha_v);
    solve();
    return Rcpp::List::create(
                        Rcpp::Named("lambda_u") = opt_lambda_u,
                        Rcpp::Named("lambda_v") = opt_lambda_v,
                        Rcpp::Named("alpha_u") = opt_alpha_u,
                        Rcpp::Named("alpha_v") = opt_alpha_v,
                        Rcpp::Named("u") = u,
                        Rcpp::Named("v") = v,
                        Rcpp::Named("d") = arma::as_scalar(u.t() * X * v));

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

arma::vec set_grid(arma::vec vec, int want_grid){
    if (want_grid == 1){ 
        return vec;
    }
    else if (want_grid == 0){
        vec.resize(1);
        vec(0) = -1;
        return vec;
    }
}

arma::vec set_bic_grid(const arma::vec &vec, int want_bic, int i){
    if(want_bic == 1){
        return vec;
    }
    else if (want_bic == 0){
        return Rcpp::NumericVector::create(vec(i));
    }
}

Rcpp::List MoMA::grid_BIC_mix(const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    int bicau,  // flags; = 0 means grid, = 01 means BIC search
    int bicav,
    int biclu,
    int biclv,
    int max_bic_iter){
    
    // If alpha_u is selected via grid search, then 
    // grid_au = alpha_u, bic_au = [-1].
    // If alpha_u is selected via nested BIC search,
    // then grid_au = [-1], bic_au = alpha_u
    arma::vec grid_lu = set_grid(lambda_u, !biclu);
    arma::vec grid_lv = set_grid(lambda_v, !biclv);
    arma::vec grid_au = set_grid(alpha_u, !bicau);
    arma::vec grid_av = set_grid(alpha_v, !bicav);

    int n_lu = grid_lu.n_elem;
    int n_lv = grid_lv.n_elem;
    int n_au = grid_au.n_elem;
    int n_av = grid_av.n_elem;

    
    Rcpp::List my_list(n_lu * n_lv * n_au * n_av);
    my_list.attr("dim") = Rcpp::NumericVector::create(
                n_au, n_lu, n_av, n_lv);

    Rcpp::List u_result;
    Rcpp::List v_result;

    for(int i = 0; i < n_au; i++){
        for(int j = 0; j < n_lu; j++){
            for(int k = 0; k < n_av; k++){
                for(int m = 0; m < n_lv; m++){

                    arma::vec bic_au = set_bic_grid(alpha_u, bicau, i);
                    arma::vec bic_lu = set_bic_grid(lambda_u, biclu, j);
                    arma::vec bic_av = set_bic_grid(alpha_v, bicav, k);
                    arma::vec bic_lv = set_bic_grid(lambda_v, biclv, m);

                    tol = 1;
                    iter = 0;
                    arma::vec oldu;
                    arma::vec oldv;
                    while(tol > EPS && iter < MAX_ITER && iter < max_bic_iter){
                        iter++;
                        oldu = u;
                        oldv = v;

                        // choose lambda/alpha_u
                        MoMALogger::debug("Start u search.");
                        u_result = bicsr_u.search(X*v,u,bic_au,bic_lu);
                        u = Rcpp::as<Rcpp::NumericVector>(u_result["vector"]);

                        MoMALogger::debug("Start v search.");
                        v_result = bicsr_v.search(X.t()*u,v,bic_av,bic_lv);
                        v = Rcpp::as<Rcpp::NumericVector>(v_result["vector"]);

                        tol = norm(oldu - u) / norm(oldu) + norm(oldv - v) / norm(oldv);
                        MoMALogger::message("Finish BIC search outer loop. (iter, tol) = (") 
                                    << iter << "," << tol << "), "
                                    << "(bic_u, bic_v) = (" 
                                    << (double)u_result["bic"] << "," 
                                    << (double)v_result["bic"] << ")";
                    }
                   
                    my_list(n_lu * n_av * n_lv * i + 
                                   n_av * n_lv * j + 
                                          n_lv * k +
                                                 m)
                            = Rcpp::List::create(
                                Rcpp::Named("u") = u_result,
                                Rcpp::Named("v") = v_result);
                
                }
            }
        }
    }

    return my_list;
}