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

const arma::vec &set_greedy_grid(const arma::vec &grid, int want_grid){
    if (want_grid == 1){ 
        return grid;
    }
    else if (want_grid == 0){
        return MOMA_EMPTY_GRID_OF_LENGTH1;
    }
}

arma::vec set_bic_grid(const arma::vec &grid, int want_bic, int i){
    if(want_bic == 1){
        return grid;
    }
    else if (want_bic == 0){
        return grid(i) * arma::ones<arma::vec> (1);
    }
}

Rcpp::List MoMA::grid_BIC_mix(const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    int bic_search_alpha_u,  // flags; = 0 means grid, = 01 means BIC search
    int bic_search_alpha_v,
    int bic_search_lambda_u,
    int bic_search_lambda_v,
    int max_bic_iter){
    
    // If alpha_u is selected via grid search, then the variable
    // grid_au = alpha_u, bic_au_grid = [-1].
    // If alpha_u is selected via nested BIC search,
    // then grid_au = [-1], bic_au_grid = alpha_u
    const arma::vec &grid_lu = set_greedy_grid(lambda_u, !bic_search_lambda_u);
    const arma::vec &grid_lv = set_greedy_grid(lambda_v, !bic_search_lambda_v);
    const arma::vec &grid_au = set_greedy_grid(alpha_u, !bic_search_alpha_u);
    const arma::vec &grid_av = set_greedy_grid(alpha_v, !bic_search_alpha_v);

    // Test that if a grid is set to be BIC-search grid, then
    // the above code should set grid_xx to the vector [-1]
    if((bic_search_alpha_u == 1 && (grid_au.n_elem != 1 || grid_au(0) != -1))
        || (bic_search_alpha_v == 1 && (grid_av.n_elem != 1 || grid_av(0) != -1))
        || (bic_search_lambda_u == 1 && (grid_lu.n_elem != 1 || grid_lu(0) != -1))
        || (bic_search_lambda_v == 1 && (grid_lv.n_elem != 1 || grid_lv(0) != -1)) )
    {
        MoMALogger::error("Wrong grid-search grid!")
        << "grid_lu.n_elem=" << grid_lu.n_elem
        << ", grid_av.n_elem=" << grid_av.n_elem
        << ", grid_lu.n_elem" << grid_lu.n_elem
        << ", grid_lv.n_elem" << grid_lv.n_elem;
    }

    int n_lu = grid_lu.n_elem;
    int n_lv = grid_lv.n_elem;
    int n_au = grid_au.n_elem;
    int n_av = grid_av.n_elem;

    
    RcppFourDList four_d_list(n_au, n_lu, n_av, n_lv);

    // nested-BIC search returns a list that
    // contains (lambda, alpha, bic, selected vector)
    Rcpp::List u_result;
    Rcpp::List v_result;

    // to facilitate warm-start
    arma::vec oldu;
    arma::vec oldv;

    for(int i = 0; i < n_au; i++){
        for(int j = 0; j < n_lu; j++){
            for(int k = 0; k < n_av; k++){
                for(int m = 0; m < n_lv; m++){

                    arma::vec bic_au_grid = set_bic_grid(alpha_u, bic_search_alpha_u, i);
                    arma::vec bic_lu_grid = set_bic_grid(lambda_u, bic_search_lambda_u, j);
                    arma::vec bic_av_grid = set_bic_grid(alpha_v, bic_search_alpha_v, k);
                    arma::vec bic_lv_grid = set_bic_grid(lambda_v, bic_search_lambda_v, m);

                    if((bic_search_alpha_u == 0 && bic_au_grid.n_elem != 1)
                        || (bic_search_alpha_v == 0 && bic_av_grid.n_elem != 1)
                        || (bic_search_lambda_u == 0 && bic_lu_grid.n_elem != 1)
                        || (bic_search_lambda_v == 0 && bic_lv_grid.n_elem != 1) )
                    {
                            MoMALogger::error("Wrong BIC search grid!");
                    }

                    tol = 1;
                    iter = 0;

                    // We conduct 2 BIC searches over 2D grids here instead 
                    // of 4 searches over 1D grids. It's consistent with 
                    // Genevera's code
                    while(tol > EPS && iter < MAX_ITER && iter < max_bic_iter){
                        iter++;
                        oldu = u;
                        oldv = v;

                        // choose lambda/alpha_u
                        MoMALogger::debug("Start u search.");
                        u_result = bicsr_u.search(X*v, u, bic_au_grid, bic_lu_grid);
                        u = Rcpp::as<Rcpp::NumericVector>(u_result["vector"]);

                        MoMALogger::debug("Start v search.");
                        v_result = bicsr_v.search(X.t()*u, v, bic_av_grid, bic_lv_grid);
                        v = Rcpp::as<Rcpp::NumericVector>(v_result["vector"]);

                        tol = norm(oldu - u) / norm(oldu) + norm(oldv - v) / norm(oldv);
                        MoMALogger::message("Finish BIC search outer loop. (iter, tol) = (") 
                                    << iter << "," << tol << "), "
                                    << "(bic_u, bic_v) = (" 
                                    << (double)u_result["bic"] << "," 
                                    << (double)v_result["bic"] << ")";
                    }
                   
                    four_d_list.insert(Rcpp::List::create(
                                Rcpp::Named("u") = u_result,
                                Rcpp::Named("v") = v_result), i, j, k, m);
 
                }
            }
        }
    }

    return four_d_list.get_list();
}
