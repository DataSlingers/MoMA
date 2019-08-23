// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
#include "moma.h"

// Initializer for PCA
MoMA::MoMA(const arma::mat &i_X,  // Pass X_ as a reference to avoid copy
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
           double i_alpha_u,  // Smoothing levels
           double i_alpha_v,
           const arma::mat &i_Omega_u,  // Smoothing matrices
           const arma::mat &i_Omega_v,

           /*
            * Algorithm parameters:
            */
           double i_EPS,
           long i_MAX_ITER,
           double i_EPS_inner,
           long i_MAX_ITER_inner,
           std::string i_solver,
           DeflationScheme i_ds)
    : n(i_X.n_rows),
      p(i_X.n_cols),
      alpha_u(i_alpha_u),
      alpha_v(i_alpha_v),
      lambda_u(i_lambda_u),
      lambda_v(i_lambda_v),
      X(i_X),  // no copy of the data
      ds(i_ds),
      Omega_u(i_Omega_u),
      Omega_v(i_Omega_v),
      MAX_ITER(i_MAX_ITER),
      EPS(i_EPS),
      solver_u(i_solver,
               alpha_u,
               i_Omega_u,
               lambda_u,
               i_prox_arg_list_u,
               i_EPS_inner,
               i_MAX_ITER_inner,
               i_X.n_rows),
      solver_v(i_solver,
               alpha_v,
               i_Omega_v,
               lambda_v,
               i_prox_arg_list_v,
               i_EPS_inner,
               i_MAX_ITER_inner,
               i_X.n_cols)
// const reference must be passed to initializer list
{
    X_original = i_X;

    if (i_EPS >= 1 || i_EPS_inner >= 1)
    {
        MoMALogger::error("EPS or EPS_inner too large.");
    }

    bicsr_u.bind(&solver_u, &PR_solver::bic);
    bicsr_v.bind(&solver_v, &PR_solver::bic);

    MoMALogger::info("Initializing MoMA object:")
        << " lambda_u " << lambda_u << " lambda_v " << lambda_v << " alpha_u " << alpha_u
        << " alpha_v " << alpha_v << " P_u " << Rcpp::as<std::string>(i_prox_arg_list_u["P"])
        << " P_v " << Rcpp::as<std::string>(i_prox_arg_list_v["P"]) << " EPS " << i_EPS
        << " MAX_ITER " << i_MAX_ITER << " EPS_inner " << i_EPS_inner << " MAX_ITER_inner "
        << i_MAX_ITER_inner << " solver " << i_solver;
    // Step 2: Initialize to leading singular vectors
    //
    //         MoMA is a regularized SVD, which is a non-convex (bi-convex)
    //         problem, so we need to be cautious about initialization to
    //         avoid local-minima. Initialization at the SVD (global solution
    //         to the non-regularized problem) seems to be a good trade-off:
    //         for problems with little regularization, the MoMA solution will
    //         lie near the SVD solution; for problems with significant
    //         regularization the problem becomes more well-behaved and less
    //         sensitive to initialization
    initialize_uv();
    is_initialzied = true;
    is_solved      = false;  // TODO: check if alphauv == 0
};

// Initializer for LDA and CCA
MoMA::MoMA(
    // const arma::mat &X_,  // Pass X_ as a reference to avoid copy
    const arma::mat &i_X_working,
    const arma::mat &i_Y_working,
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
    double i_alpha_u,  // Smoothing levels
    double i_alpha_v,
    const arma::mat &i_Omega_u,  // Smoothing matrices
    const arma::mat &i_Omega_v,

    /*
     * Algorithm parameters:
     */
    double i_EPS,
    long i_MAX_ITER,
    double i_EPS_inner,
    long i_MAX_ITER_inner,
    std::string i_solver,
    DeflationScheme i_ds)
    : MoMA(i_X_working.t() * i_Y_working,  // the matrix on which we find pSVD
           i_lambda_u,
           i_lambda_u,
           i_prox_arg_list_u,
           i_prox_arg_list_v,
           i_alpha_u,
           i_alpha_v,
           i_Omega_u,
           i_Omega_v,
           i_EPS,
           i_MAX_ITER,
           i_EPS_inner,
           i_MAX_ITER_inner,
           i_solver,
           i_ds)
{
    if (ds == DeflationScheme::CCA)
    {
        // const matrix
        X_original = i_X_working;
        Y_original = i_Y_working;

        // deflated matrices
        X_working = i_X_working;
        Y_working = i_Y_working;
    }
    else if (ds == DeflationScheme::LDA)
    {
        MoMALogger::debug("Initializing MoMA LDA mode.");
        // const matrix
        X_original = i_X_working;
        Y_original = i_Y_working;

        // deflated matrices
        X_working = i_X_working;
        // we do not need Y_working
        // since Y is an indicator matrix
        MoMALogger::debug(" X_target = \n") << X;
    }
    else
    {
        MoMALogger::error("Error in MoMA LDA initialzation.");
    }
};

arma::vec normalize(const arma::vec &u)
{
    arma::vec res = u;
    double mn     = arma::norm(u);
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

int MoMA::deflate()
{
    if (is_solved != true)
    {
        MoMALogger::error("Please call `MoMA::solve` before `MoMA::deflate`.");
    }
    if (ds == DeflationScheme::PCA_Hotelling)
    {
        double d = arma::as_scalar(u.t() * X * v);
        MoMALogger::debug("Deflating:\n")
            << "\nX = \n"
            << X << "u^T = " << u.t() << "v^T = " << v.t() << "d = u^TXv = " << d;

        if (d <= 0.0)
        {
            MoMALogger::error("Cannot deflate by non-positive factor.");
        }
        X = X - d * u * v.t();
        // Re-initialize u and v after deflation
        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::PCA_Schur_Complement)
    {
        double d = arma::as_scalar(u.t() * X * v);
        if (d <= 0.0)
        {
            MoMALogger::error("Error in Schur complement: devision by zero.");
        }

        // No need to scale u and v
        X = X - (X * v) * (u.t() * X) / d;

        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::PCA_Projection)
    {
        arma::mat eye_u(n, n, arma::fill::eye);
        arma::mat eye_v(p, p, arma::fill::eye);

        arma::vec u_unit = normalize(u);
        arma::vec v_unit = normalize(v);

        X = (eye_u - u_unit * u_unit.t()) * X * (eye_v - v_unit * v_unit.t());

        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::CCA)
    {
        double u_norm = arma::norm(u);
        double v_norm = arma::norm(v);

        if (u_norm == 0.0 || v_norm == 0.0)
        {
            MoMALogger::error("Zero singular vecters in MoMA::deflate.");
        }

        // u and v are scores
        // "cv" = canonical variates
        arma::mat X_cv   = X_working * u;
        arma::mat Y_cv   = Y_working * v;
        double norm_X_cv = arma::norm(X_cv);
        double norm_Y_cv = arma::norm(Y_cv);

        // subtract cv's out of X_working and Y_working
        X_working = X_working - 1 / (norm_X_cv * norm_X_cv) * X_cv * X_cv.t() * X_working;
        Y_working = Y_working - 1 / (norm_Y_cv * norm_Y_cv) * Y_cv * Y_cv.t() * Y_working;

        X = X_working.t() * Y_working;

        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::LDA)
    {
        double u_norm = arma::norm(u);
        double v_norm = arma::norm(v);

        if (u_norm == 0.0 || v_norm == 0.0)
        {
            MoMALogger::error("Zero singular vecters in MoMA::deflate.");
        }

        // u and v are scores
        // "cv" = canonical variates
        arma::mat X_cv   = X_working * u;
        double norm_X_cv = arma::norm(X_cv);

        // subtract cv's out of X_working only
        X_working = X_working - 1 / (norm_X_cv * norm_X_cv) * X_cv * X_cv.t() * X_working;

        X = X_working.t() * Y_original;

        initialize_uv();
        return 0;
    }
    else
    {
        MoMALogger::error("Wrong defaltion scheme.");
    }
}

// Dependence on MoMA's internal states: MoMA::X, MoMA::u, MoMA::v, MoMA::alpha_u/v,
// MoMA::lambda_u/v
// After calling MoMA::solve(), MoMA::u and MoMA::v become the solution to the penalized regression.
void MoMA::solve()
{
    double tol = 1;
    int iter   = 0;
    arma::vec oldu;
    arma::vec oldv;
    while (tol > EPS && iter < MAX_ITER)
    {
        iter++;
        oldu = u;
        oldv = v;

        u = solver_u.solve(X * v, u);
        v = solver_v.solve(X.t() * u, v);

        double scale_u = arma::norm(oldu) == 0.0 ? 1 : arma::norm(oldu);
        double scale_v = arma::norm(oldv) == 0.0 ? 1 : arma::norm(oldv);

        tol = arma::norm(oldu - u) / scale_u + arma::norm(oldv - v) / scale_v;
        MoMALogger::debug("Real-time PG loop info:  (iter, tol) = (") << iter << ", " << tol << ")";
    }

    MoMALogger::info("Finish PG loop. Total iter = ") << iter;
    check_convergence(iter, tol);
    is_solved = true;
}

double MoMA::evaluate_loss()
{
    if (!is_solved)
    {
        MoMALogger::error("Please call MoMA::solve first before MoMA::evaluate_loss.");
    }
    double u_ellipsoid_constraint = arma::as_scalar(u.t() * u + alpha_u * u.t() * Omega_u * u);
    double v_ellipsoid_constraint = arma::as_scalar(v.t() * v + alpha_v * v.t() * Omega_v * v);

    if ((std::abs(u_ellipsoid_constraint) > MOMA_FLOATPOINT_EPS &&
         std::abs(u_ellipsoid_constraint - 1.0) > MOMA_FLOATPOINT_EPS) ||
        (std::abs(v_ellipsoid_constraint) > MOMA_FLOATPOINT_EPS &&
         std::abs(v_ellipsoid_constraint - 1.0) > MOMA_FLOATPOINT_EPS))
    {
        MoMALogger::error("Ellipse constraint is not met.");
    }

    MoMALogger::error("MoMA::evaluate_loss it not implemented yet.");
    return 0;  // TODO: Implement it.
}

int MoMA::initialize_uv()
{
    // TODO: we can set MoMA::u and MoMA::v to
    // the solution of pSVD with only smoothness constraints.

    // Set MoMA::v, MoMA::u as leading SVs of X
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    v              = V.col(0);
    u              = U.col(0);
    is_initialzied = true;
    return 0;
}

int MoMA::check_convergence(int iter, double tol)
{
    if (iter >= MAX_ITER || tol > EPS)
    {
        MoMALogger::warning("No convergence in MoMA!")
            << " lambda_u " << lambda_u << " lambda_v " << lambda_v << " alpha_u " << alpha_u
            << " alpha_v " << alpha_v;
    }
    return 0;
}

// Note it does not change MoMA::u and MoMA::v
int MoMA::set_penalty(double newlambda_u, double newlambda_v, double newalpha_u, double newalpha_v)
{
    solver_u.set_penalty(newlambda_u, newalpha_u);
    solver_v.set_penalty(newlambda_v, newalpha_v);

    // NOTE: We must keep the alpha's and lambda's up-to-date
    // in both MoMA and solve_u/v
    if (std::abs(alpha_u - newalpha_u) > MOMA_FLOATPOINT_EPS ||
        std::abs(alpha_v - newalpha_v) > MOMA_FLOATPOINT_EPS ||
        std::abs(lambda_v - newlambda_v) > MOMA_FLOATPOINT_EPS ||
        std::abs(lambda_u - newlambda_u) > MOMA_FLOATPOINT_EPS)
    {
        // update internal states of MoMA
        is_solved = false;
        alpha_u   = newalpha_u;
        alpha_v   = newalpha_v;
        lambda_u  = newlambda_u;
        lambda_v  = newlambda_v;
    }

    return 0;
}

int MoMA::reset_X()
{
    if (ds == DeflationScheme::PCA_Hotelling || ds == DeflationScheme::PCA_Schur_Complement ||
        ds == DeflationScheme::PCA_Projection)
    {
        X = X_original;
        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::CCA)
    {
        X_working = X_original;
        Y_working = Y_original;
        X         = X_working.t() * Y_working;
        initialize_uv();
        return 0;
    }
    else if (ds == DeflationScheme::LDA)
    {
        X_working = X_original;
        X         = X_working.t() * Y_original;
    }
    else
    {
        MoMALogger::error("MoMA::reset_X for other modes not implemented.");
    }
}
