// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma_solver.h"

// A handle class
// Penalized regression solver
// min_u || y - u || + lambda * P(u) s.t. || u ||_S <= 1
// S = I + alpha * Omega
arma::vec _PR_solver::normalize(const arma::vec &u){
    arma::vec res = u;
    double mn = mat_norm(u,S,I);
    if(mn > 0){
        res /= mn;
    }else{
        res.zeros();
    }
    return res;
}

_PR_solver::_PR_solver(
        double i_alpha, const arma::mat &i_Omega, double i_lambda,
        const std::string &sparsity_string, double gamma,
        const arma::vec &group, const arma::mat &w,
        bool ADMM, bool acc, double prox_eps, bool nonneg,
        double i_EPS, int i_MAX_ITER):
        lambda(i_lambda),
        alpha(i_alpha),
        Omega(i_Omega),
        p(sparsity_string,gamma,group,w,ADMM,acc,prox_eps,nonneg),
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
    arma::vec old;    // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        old = u;

        u = g(u,y,grad_step_size,S,I);
        u = p(u,prox_step_size);

        tol = arma::norm(u - old) / arma::norm(old);
        MoMALogger::debug("No.") << iter << "--"<< "%change " << tol;
    }
    u = normalize(u);
    
    return u;
}

arma::vec FISTA::solve(arma::vec y, const arma::vec &start_point){
    if(start_point.n_elem != S.n_cols || y.n_elem != S.n_cols){
        MoMALogger::error("Wroing dimension in PRsolver::solve");
    }
    tol = 1;
    iter = 0;
    arma::vec u = start_point;
    arma::vec old;    // store working result

    double t = 1;   
    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        old = u;
        double oldt = t;
        t = 0.5 * (1 + std::sqrt(1 + 4 * oldt*oldt));

        u = g(u,y,grad_step_size,S,I);
        u = p(u,prox_step_size);
        u = u + (oldt - 1) / t * (u - old);

        tol = arma::norm(u - old) / arma::norm(old);
        MoMALogger::debug("No.") << iter << "--"<< "%change " << tol;
    }
    u = normalize(u);
    
    return u;
}

arma::vec OneStepISTA::solve(arma::vec y, const arma::vec &start_point){
    if(start_point.n_elem != S.n_cols || y.n_elem != S.n_cols){
        MoMALogger::error("Wroing dimension in PRsolver::solve");
    }
    tol = 1;
    iter = 0;
    arma::vec u = start_point;
    arma::vec old;    // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        old = u;

        u = g(u,y,grad_step_size,S,I);
        u = p(u,prox_step_size);
        u = normalize(u);

        tol = arma::norm(u - old) / arma::norm(old);
        MoMALogger::debug("No.") << iter << "--"<< "%change " << tol;
    }
    
    return u;
}

PR_solver::PR_solver(
    const std::string &algorithm_string,
    double i_alpha,const arma::mat &i_Omega,
    double i_lambda,
    const std::string &sparsity_string, double gamma,
    const arma::vec &group,
    const arma::mat &w, bool ADMM, bool acc, double prox_eps,
    bool nonneg,
    double i_EPS,int i_MAX_ITER){

    if (algorithm_string.compare("ISTA") == 0){
        prs = new ISTA(
                    i_alpha,i_Omega,i_lambda,sparsity_string,
                    gamma,group,w,ADMM,acc,prox_eps,nonneg,
                    i_EPS,i_MAX_ITER);
    }
    else if (algorithm_string.compare("FISTA") == 0){
        prs =  new FISTA(
                    i_alpha,i_Omega,i_lambda,sparsity_string,
                    gamma,group,w,ADMM,acc,prox_eps,nonneg,
                    i_EPS,i_MAX_ITER);
    }
    else if (algorithm_string.compare("FISTA") == 0){
        prs =  new OneStepISTA(
                    i_alpha,i_Omega,i_lambda,sparsity_string,
                    gamma,group,w,ADMM,acc,prox_eps,nonneg,
                    i_EPS,i_MAX_ITER);
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