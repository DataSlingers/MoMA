// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma_solver.h"

// A handle class
// Penalized regression solver
// min_u || y - u || + lambda * P(u) s.t. || u ||_S <= 1
// S = I + alpha * Omega
arma::vec _PR_solver::normalize(const arma::vec &u){
    arma::vec res = u;
    double mn = I? arma::norm(u) : arma::as_scalar(arma::sqrt(u.t() * S * u));
    if(mn > 0){
        res /= mn;
    }else{
        res.zeros();
    }
    return res;
}

int _PR_solver::check_cnvrg(){
    if(iter >= MAX_ITER){
        MoMALogger::warning("No convergence in _PR_solver!");
    }
    return 0;
}

_PR_solver::_PR_solver(
        double i_alpha, const arma::mat &i_Omega,
        double i_lambda, Rcpp::List prox_arg_list,
        double i_EPS, int i_MAX_ITER, int i_dim):
        dim(i_dim),
        lambda(i_lambda),
        alpha(i_alpha),
        Omega(i_Omega),         // reference to the matrix on the R side, no extra copy
        p(prox_arg_list,i_dim),
        EPS(i_EPS),
        MAX_ITER(i_MAX_ITER){

    // Step 1b: Calculate leading eigenvalues of smoothing matrices
    //          -> used for prox gradient step sizes
    S.eye(arma::size(Omega));
    S += alpha * Omega;
    L = arma::eig_sym(S).max() + MOMA_EIGENVALUE_REGULARIZATION;

    grad_step_size = 1 / L;
    prox_step_size = lambda / L;
    I = (alpha == 0.0);
    tol = 1;
    iter = 0;
}

arma::vec _PR_solver::g(
                    const arma::vec &v, const arma::vec &y,
                    double step_size, const arma::mat &S, bool I){
    arma::vec res;
    if(I){
        res = v + step_size * (y - v);
    }else{
        res = v + step_size * (y - S*v);
    }
    return res;
}

int _PR_solver::reset(double new_lambda, double new_alpha){
    if(new_alpha != alpha){
        // avoid re-calculating L
        // enter only when alpha is changed
        S.eye(arma::size(S));
        S += new_alpha * Omega;
        L = arma::eig_sym(S).max() + MOMA_EIGENVALUE_REGULARIZATION;

        grad_step_size = 1 / L;
        prox_step_size = new_lambda / L;
    }
    else if(lambda != new_lambda){
        prox_step_size = new_lambda / L;
    }
    lambda = new_lambda;
    alpha = new_alpha;
    I = (new_alpha == 0.0);
    return 0;
}

arma::vec ISTA::solve(arma::vec y, const arma::vec &start_point){
    if(start_point.n_elem != S.n_cols || y.n_elem != S.n_cols){
        MoMALogger::error("Wroing dimension in PRsolver::solve");
    }
    tol = 1;
    iter = 0;
    arma::vec u = start_point;
    arma::vec oldu;    // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;

        u = g(u,y,grad_step_size,S,I);
        u = p(u,prox_step_size);

        tol = arma::norm(u - oldu) / arma::norm(oldu);
        if(iter % 1000 == 0){
            MoMALogger::debug("No.") << iter << "--"<< tol;
        }
    }
    u = normalize(u);
    
    MoMALogger::debug("Inner loop No.") << iter << "--" << tol;
    check_cnvrg();
    return u;
}

arma::vec FISTA::solve(arma::vec y, const arma::vec &start_point){
    if(start_point.n_elem != S.n_cols || y.n_elem != S.n_cols){
        MoMALogger::error("Wroing dimension in PRsolver::solve");
    }
    tol = 1;
    iter = 0;
    arma::vec u = start_point;
    arma::vec newu = start_point;
    arma::vec oldu;    // store working result

    double t = 1;   
    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;
        double oldt = t;
        t = 0.5 * (1 + std::sqrt(1 + 4 * oldt*oldt));

        u = g(newu,y,grad_step_size,S,I);
        u = p(u,prox_step_size);
        newu = u + (oldt - 1) / t * (u - oldu);

        tol = arma::norm(u - oldu) / arma::norm(oldu);
        if(iter % 1000 == 0){
            MoMALogger::debug("No.") << iter << "--"<< tol;
        }
    }
    u = normalize(u);
    
    check_cnvrg();
    MoMALogger::debug("Inner loop No.") << iter << "--" << tol;
    return u;
}

arma::vec OneStepISTA::solve(arma::vec y, const arma::vec &start_point){
    if(start_point.n_elem != S.n_cols || y.n_elem != S.n_cols){
        MoMALogger::error("Wroing dimension in PRsolver::solve");
    }
    tol = 1;
    iter = 0;
    arma::vec u = start_point;
    arma::vec oldu;    // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;

        u = g(u,y,grad_step_size,S,I);
        u = p(u,prox_step_size);
        u = normalize(u);

        tol = arma::norm(u - oldu) / arma::norm(oldu);
        if(iter % 1000 == 0){
            MoMALogger::debug("No.") << iter << "--"<< tol;
        }
    }
    
    check_cnvrg();
    MoMALogger::debug("Inner loop No.") << iter << "--" << tol;
    return u;
}

PR_solver::PR_solver(
    const std::string &algorithm_string,
    double i_alpha, const arma::mat &i_Omega,
    double i_lambda, Rcpp::List prox_arg_list,
    double i_EPS, int i_MAX_ITER, int dim){

    if (algorithm_string.compare("ISTA") == 0){
        prs = new ISTA(
                    i_alpha,i_Omega,
                    i_lambda,prox_arg_list,
                    i_EPS,i_MAX_ITER,dim);
    }
    else if (algorithm_string.compare("FISTA") == 0){
        prs =  new FISTA(
                    i_alpha,i_Omega,
                    i_lambda,prox_arg_list,
                    i_EPS,i_MAX_ITER,dim);
    }
    else if (algorithm_string.compare("ONESTEPISTA") == 0){
        prs =  new OneStepISTA(
                    i_alpha,i_Omega,
                    i_lambda,prox_arg_list,
                    i_EPS,i_MAX_ITER,dim);
    }
    else{
        MoMALogger::error("Your choice of algorithm not provided: ") << algorithm_string;
    }
};

arma::vec PR_solver::solve(arma::vec y, const arma::vec &start_point){
    return (*prs).solve(y,start_point);
}
    
int PR_solver::reset(double new_lambda, double new_alpha){
    return (*prs).reset(new_lambda,new_alpha);
}
