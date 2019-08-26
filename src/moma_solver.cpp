// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
#include "moma_solver.h"

// A handle class
// Penalized regression solver
// min_u || y - u || + lambda * P(u) s.t. || u ||_S <= 1
// S = I + alpha * Omega
arma::vec _PR_solver::normalize(const arma::vec &u)
{
    arma::vec res = u;
    double mn     = is_S_idmat ? arma::norm(u) : arma::as_scalar(arma::sqrt(u.t() * S * u));
    if (mn > 0)
    {
        res /= mn;
    }
    else
    {
        res.zeros();
    }
    return res;
}

void _PR_solver::check_convergence(int iter, double tol)
{
    if (iter >= MAX_ITER || tol > EPS)
    {
        MoMALogger::warning("No convergence in _PR_solver!");
    }
}

_PR_solver::_PR_solver(double i_alpha,
                       const arma::mat &i_Omega,
                       double i_lambda,
                       Rcpp::List prox_arg_list,
                       double i_EPS,
                       int i_MAX_ITER,
                       int i_dim)
    : dim(i_dim),
      lambda(i_lambda),
      alpha(i_alpha),
      Omega(i_Omega),  // reference to the matrix on the R side, no extra copy
      p(prox_arg_list, i_dim),
      EPS(i_EPS),
      MAX_ITER(i_MAX_ITER)
{
    // Step 1b: Calculate leading eigenvalues of smoothing matrices
    //          -> used for prox gradient step sizes
    S.eye(arma::size(Omega));
    S += alpha * Omega;
    L = arma::eig_sym(S).max() + MOMA_EIGENVALUE_REGULARIZATION;

    grad_step_size = 1 / L;
    prox_step_size = lambda / L;
    is_S_idmat     = (alpha == 0.0);
}

arma::vec _PR_solver::g(const arma::vec &v,
                        const arma::vec &y,
                        double step_size,
                        const arma::mat &S,
                        bool is_S_idmat)
{
    arma::vec res;
    if (is_S_idmat)
    {
        res = v + step_size * (y - v);
    }
    else
    {
        res = v + step_size * (y - S * v);
    }
    return res;
}

int _PR_solver::set_penalty(double new_lambda, double new_alpha)
{
    if (new_alpha != alpha)
    {
        // avoid re-calculating L
        // enter only when alpha is changed
        S.eye(arma::size(S));
        S += new_alpha * Omega;
        L = arma::eig_sym(S).max() + MOMA_EIGENVALUE_REGULARIZATION;

        grad_step_size = 1 / L;
        prox_step_size = new_lambda / L;
    }
    else if (lambda != new_lambda)
    {
        prox_step_size = new_lambda / L;
    }
    lambda     = new_lambda;
    alpha      = new_alpha;
    is_S_idmat = (new_alpha == 0.0);
    return 0;
}

double _PR_solver::bic(arma::vec y, const arma::vec &est)
{
    // Find out the bic of the following estimator:
    // argmin_x || y - x || + lambda P(x) s.t. || x ||_S <= 1, S = I + alpah *
    // Omega. It is approximated by log{1/2n * squared error} + log(n)/n * df +
    // const NOTE: We ignore the effect of smoothing matrix Omega.

    // Ref: Proposition 2 in
    // Allen, Genevera I.
    // "Sparse and functional principal components analysis."
    // arXiv preprint arXiv:1309.2895 (2013).
    double res = arma::norm(y - est);
    double df  = p.df(est);
    MoMALogger::debug("(RES, DF) = （") << res << ", " << df << ").";
    if (res == 0.0)
    {
        MoMALogger::warning("BIC = -infty due to zero resdiual.");
    }
    return std::log(res * res / dim) + std::log(dim) / dim * df;  // ignore some constants here
}

arma::vec ISTA::solve(arma::vec y, const arma::vec &start_point)
{
    if (start_point.n_elem != S.n_cols || y.n_elem != S.n_cols)
    {
        MoMALogger::error("Wrong dimension in PRsolver::solve:")
            << start_point.n_elem << ":" << S.n_cols;
    }
    double tol  = 1;
    int iter    = 0;
    arma::vec u = start_point;
    arma::vec oldu;  // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;

        u = g(u, y, grad_step_size, S, is_S_idmat);
        u = p(u, prox_step_size);

        double old_norm  = arma::norm(oldu);
        double diff_norm = arma::norm(u - oldu);

        if (old_norm != 0.0)
            tol = diff_norm / old_norm;
        else
        {
            tol = diff_norm;
        }
        if (iter % 1000 == 0)
        {
            MoMALogger::debug("Solving PR: No.") << iter << "--" << tol;
        }
    }
    u = normalize(u);

    MoMALogger::debug("Finish solving PR: (total_iter, tol) = ")
        << "(" << iter << "," << tol << ")";
    check_convergence(iter, tol);
    return u;
}

arma::vec FISTA::solve(arma::vec y, const arma::vec &start_point)
{
    if (start_point.n_elem != S.n_cols || y.n_elem != S.n_cols)
    {
        MoMALogger::error("Wrong dimension in PRsolver::solve");
    }
    double tol     = 1;
    int iter       = 0;
    arma::vec u    = start_point;
    arma::vec newu = start_point;
    arma::vec oldu;  // store working result

    double t = 1;
    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu        = u;
        double oldt = t;
        t           = 0.5 * (1 + std::sqrt(1 + 4 * oldt * oldt));

        u    = g(newu, y, grad_step_size, S, is_S_idmat);
        u    = p(u, prox_step_size);
        newu = u + (oldt - 1) / t * (u - oldu);

        double old_norm  = arma::norm(oldu);
        double diff_norm = arma::norm(u - oldu);
        if (old_norm != 0.0)
            tol = diff_norm / old_norm;
        else
        {
            tol = diff_norm;
        }
        if (iter % 1000 == 0)
        {
            MoMALogger::debug("Solving PR: No.") << iter << "--" << tol;
        }
    }
    u = normalize(u);

    check_convergence(iter, tol);
    MoMALogger::debug("Finish solving PR: (total_iter, tol) = ")
        << "(" << iter << "," << tol << ")";
    return u;
}

arma::vec OneStepISTA::solve(arma::vec y, const arma::vec &start_point)
{
    if (start_point.n_elem != S.n_cols || y.n_elem != S.n_cols)
    {
        MoMALogger::error("Wrong dimension in PRsolver::solve");
    }
    double tol  = 1;
    int iter    = 0;
    arma::vec u = start_point;
    arma::vec oldu;  // store working result

    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;

        u = g(u, y, grad_step_size, S, is_S_idmat);
        u = p(u, prox_step_size);
        u = normalize(u);

        double old_norm  = arma::norm(oldu);
        double diff_norm = arma::norm(u - oldu);
        if (old_norm != 0.0)
            tol = diff_norm / old_norm;
        else
        {
            tol = diff_norm;
        }

        if (iter % 1000 == 0)
        {
            MoMALogger::debug("Solving PR: No.") << iter << "--" << tol;
        }
    }

    check_convergence(iter, tol);
    MoMALogger::debug("Finish solving PR: (total_iter, tol) = ")
        << "(" << iter << "," << tol << ")";
    return u;
}

PR_solver::PR_solver(const std::string &algorithm_string,
                     double i_alpha,
                     const arma::mat &i_Omega,
                     double i_lambda,
                     Rcpp::List prox_arg_list,
                     double i_EPS,
                     int i_MAX_ITER,
                     int dim)
{
    if (algorithm_string.compare("ISTA") == 0)
    {
        prs = new ISTA(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim);
    }
    else if (algorithm_string.compare("FISTA") == 0)
    {
        prs = new FISTA(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim);
    }
    else if (algorithm_string.compare("ONESTEPISTA") == 0)
    {
        prs = new OneStepISTA(i_alpha, i_Omega, i_lambda, prox_arg_list, i_EPS, i_MAX_ITER, dim);
    }
    else
    {
        MoMALogger::error("Your choice of algorithm not provided: ") << algorithm_string;
    }
};

arma::vec PR_solver::solve(arma::vec y, const arma::vec &start_point)
{
    return (*prs).solve(y, start_point);
}

int PR_solver::set_penalty(double new_lambda, double new_alpha)
{
    return (*prs).set_penalty(new_lambda, new_alpha);
}

double PR_solver::bic(arma::vec y, const arma::vec &est)
{
    return (*prs).bic(y, est);
}
