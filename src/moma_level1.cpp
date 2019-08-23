// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
#include "moma.h"

// auxiliary functions for MoMA::grid_BIC_mix
const arma::vec &construct_grid_for_search(const arma::vec &grid, SelectionScheme ss)
{
    if (ss == SelectionScheme::grid)
    {
        return grid;
    }
    else if (ss == SelectionScheme::BIC)
    {
        return MOMA_EMPTY_GRID_OF_LENGTH1;
    }
}

arma::vec construct_grid_no_search(const arma::vec &grid, SelectionScheme ss, int i)
{
    if (ss == SelectionScheme::BIC)
    {
        return grid;
    }
    else if (ss == SelectionScheme::grid)
    {
        return grid(i) * arma::ones<arma::vec>(1);
    }
}

// 1. Return a list of two lists
// u_result
// -- u_result["lambda"] = opt_lambda_u,
// -- u_result["alpha"] = opt_alpha_u,
// -- u_result["vector"] = working_selected_u,
// -- u_result["bic"] = minbic_u,
// v_result = same as u_result
// 2. Dependence on MoMA's internal states: MoMA::X.
// 3. After calling MoMA::criterion_search, if final_run = true, then MoMA::u and MoMA::v become the
// solution evaluated at the chosen penalty, using leading SVs as start points, and MoMA::alpha_u/v,
// MoMA::lambda_u/v become the chosen penalty. If final_run = false, then MoMA::u, MoMA::v,
// MoMA::alpha_u/v and MoMA::lambda_u/v remain unchanged. MoMA::X always remains unchanged.
Rcpp::List MoMA::criterion_search(const arma::vec &bic_au_grid,
                                  const arma::vec &bic_lu_grid,
                                  const arma::vec &bic_av_grid,
                                  const arma::vec &bic_lv_grid,
                                  arma::vec initial_u,
                                  arma::vec initial_v,
                                  double EPS_bic,
                                  int max_bic_iter,
                                  bool final_run)
{
    double tol = 1;
    int iter   = 0;

    int n_au = bic_au_grid.n_elem;
    int n_av = bic_av_grid.n_elem;
    int n_lu = bic_lu_grid.n_elem;
    int n_lv = bic_lv_grid.n_elem;

    // u_/v_result is a list of the following contents
    // Rcpp::Named("lambda") = opt_lambda_u, the chosen lambda
    // Rcpp::Named("alpha") = opt_alpha_u, the chosen alpha
    // Rcpp::Named("vector") = working_selected_u,
    // Rcpp::Named("bic") = minbic_u,
    Rcpp::List u_result;
    Rcpp::List v_result;

    // to check convergence of nested-BIC
    arma::vec oldu;
    arma::vec oldv;

    // Usually initial_u/_v are leading SVs of X.
    // See MoMA::grid_BIC_mix.
    arma::vec curu = initial_u;
    arma::vec curv = initial_v;

    // We conduct 2 BIC searches over 2D grids here instead
    // of 4 searches over 1D grids. It's consistent with
    // Genevera's code.
    if (n_au > 1 || n_av > 1 || n_lu > 1 || n_lv > 1)
    {
        while (tol > EPS_bic && iter < max_bic_iter)
        {
            iter++;
            oldu = curu;
            oldv = curv;

            // choose lambda/alpha_u
            MoMALogger::debug("Start u search.");
            u_result = bicsr_u.search(X * curv, curu, bic_au_grid, bic_lu_grid);
            curu     = Rcpp::as<Rcpp::NumericVector>(u_result["vector"]);

            MoMALogger::debug("Start v search.");
            v_result = bicsr_v.search(X.t() * curu, curv, bic_av_grid, bic_lv_grid);
            curv     = Rcpp::as<Rcpp::NumericVector>(v_result["vector"]);

            double scale_u = arma::norm(oldu) == 0.0 ? 1 : arma::norm(oldu);
            double scale_v = arma::norm(oldv) == 0.0 ? 1 : arma::norm(oldv);

            tol = arma::norm(oldu - u) / scale_u + arma::norm(oldv - v) / scale_v;
            MoMALogger::debug("Finish nested greedy BIC search outer loop. (iter, tol) = (")
                << iter << "," << tol << "), "
                << "(bic_u, bic_v) = (" << (double)u_result["bic"] << "," << (double)v_result["bic"]
                << ")";
        }
    }
    else
    {
        u_result = Rcpp::List::create(
            Rcpp::Named("lambda") = bic_lu_grid(0), Rcpp::Named("alpha") = bic_au_grid(0),
            Rcpp::Named("vector") = initial_u, Rcpp::Named("bic") = -MOMA_INFTY);
        v_result = Rcpp::List::create(
            Rcpp::Named("lambda") = bic_lv_grid(0), Rcpp::Named("alpha") = bic_av_grid(0),
            Rcpp::Named("vector") = initial_v, Rcpp::Named("bic") = -MOMA_INFTY);
        MoMALogger::debug("Deprecated BIC grid. Skip searching.");
    }

    // A final run on the selected parameter
    if (final_run)
    {
        double opt_lambda_u = u_result["lambda"];
        double opt_lambda_v = v_result["lambda"];
        double opt_alpha_u  = u_result["alpha"];
        double opt_alpha_v  = v_result["alpha"];
        MoMALogger::message("Start a final run on the chosen parameters.")
            << "[av, au, lu, lv] = [" << opt_alpha_v << ", " << opt_alpha_u << ", " << opt_lambda_u
            << ", " << opt_lambda_v << "]";

        set_penalty(opt_lambda_u, opt_lambda_v, opt_alpha_u, opt_alpha_v);
        initialize_uv();
        solve();  // Use MoMA::u and MoMA::v as start points

        u_result["vector"] = u;
        v_result["vector"] = v;

        // NOTE: we do not update bic for the new u and v.
    }

    return Rcpp::List::create(Rcpp::Named("u_result") = u_result,
                              Rcpp::Named("v_result") = v_result);
}

// Return a 5-D list. Each element is a list of
// following format:
// Rcpp::Named("u") = u_result, which is a list with
//     Rcpp::Named("lambda") = opt_lambda_u,
//     Rcpp::Named("alpha") = opt_alpha_u,
//     Rcpp::Named("vector") = working_selected_u,
//     Rcpp::Named("bic") = minbic_u
// Rcpp::Named("v") = v_result, same as "u"
// Rcpp::Named("k") = pc, the pc-th SVs
// Rcpp::Named("X"), the matrix with which we solve for u and v
Rcpp::List MoMA::grid_BIC_mix(const arma::vec &alpha_u,
                              const arma::vec &alpha_v,
                              const arma::vec &lambda_u,
                              const arma::vec &lambda_v,
                              SelectionScheme select_scheme_alpha_u,
                              SelectionScheme select_scheme_alpha_v,
                              SelectionScheme select_scheme_lambda_u,
                              SelectionScheme select_scheme_lambda_v,
                              int max_bic_iter,
                              int rank)
{

    if (rank <= 0)
    {
        MoMALogger::error("rank in MoMA::grid_BIC_mix should >= 1.");
    }

    // If alpha_u is selected via grid search, then the variable
    // grid_au = alpha_u, bic_au_grid = [-1].
    // If alpha_u is selected via nested BIC search,
    // then grid_au = [-1], bic_au_grid = alpha_u
    const arma::vec &grid_lu = construct_grid_for_search(lambda_u, select_scheme_lambda_u);
    const arma::vec &grid_lv = construct_grid_for_search(lambda_v, select_scheme_lambda_v);
    const arma::vec &grid_au = construct_grid_for_search(alpha_u, select_scheme_alpha_u);
    const arma::vec &grid_av = construct_grid_for_search(alpha_v, select_scheme_alpha_v);

    int n_lambda_u = grid_lu.n_elem;
    int n_lambda_v = grid_lv.n_elem;
    int n_alpha_u  = grid_au.n_elem;
    int n_alpha_v  = grid_av.n_elem;

    RcppFiveDList five_d_list(n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v, rank);

    for (int i = 0; i < n_alpha_u; i++)
    {
        for (int j = 0; j < n_lambda_u; j++)
        {
            for (int k = 0; k < n_alpha_v; k++)
            {
                for (int m = 0; m < n_lambda_v; m++)
                {
                    arma::vec bic_au_grid =
                        construct_grid_no_search(alpha_u, select_scheme_alpha_u, i);
                    arma::vec bic_lu_grid =
                        construct_grid_no_search(lambda_u, select_scheme_lambda_u, j);
                    arma::vec bic_av_grid =
                        construct_grid_no_search(alpha_v, select_scheme_alpha_v, k);
                    arma::vec bic_lv_grid =
                        construct_grid_no_search(lambda_v, select_scheme_lambda_v, m);

                    reset_X();
                    for (int pc = 0; pc < rank; pc++)
                    {
                        // MoMA::criterion_search returns a list containing two lists:
                        // u_result =
                        // Rcpp::Named("lambda") = opt_lambda_u,
                        // Rcpp::Named("alpha") = opt_alpha_u,
                        // Rcpp::Named("vector") = working_selected_u,
                        // Rcpp::Named("bic") = minbic_u.
                        // v_result = ... (same as above)
                        Rcpp::List result = criterion_search(
                            bic_au_grid, bic_lu_grid, bic_av_grid, bic_lv_grid, u, v, EPS,
                            max_bic_iter);  // u, v are leading SVs of MoMA::X
                        Rcpp::List u_result = result["u_result"];
                        Rcpp::List v_result = result["v_result"];

                        arma::vec curu = Rcpp::as<Rcpp::NumericVector>(u_result["vector"]);
                        arma::vec curv = Rcpp::as<Rcpp::NumericVector>(v_result["vector"]);
                        double d       = arma::as_scalar(curu.t() * X * curv);

                        Rcpp::List wrap_up;
                        if (ds == DeflationScheme::PCA_Hotelling ||
                            ds == DeflationScheme::PCA_Schur_Complement ||
                            ds == DeflationScheme::PCA_Projection)
                        {
                            wrap_up = Rcpp::List::create(
                                Rcpp::Named("u") = u_result, Rcpp::Named("v") = v_result,
                                Rcpp::Named("k") = pc, Rcpp::Named("X") = X, Rcpp::Named("d") = d);
                        }
                        else if (ds == DeflationScheme::CCA)
                        {
                            wrap_up = Rcpp::List::create(
                                Rcpp::Named("u") = u_result, Rcpp::Named("v") = v_result,
                                Rcpp::Named("k") = pc, Rcpp::Named("X") = X_working,
                                Rcpp::Named("Y") = Y_working,  // an extra "Y" element
                                Rcpp::Named("d") = d);
                        }
                        else if (ds == DeflationScheme::LDA)
                        {
                            wrap_up = Rcpp::List::create(
                                Rcpp::Named("u") = u_result, Rcpp::Named("v") = v_result,
                                Rcpp::Named("k") = pc, Rcpp::Named("X") = X_working,
                                Rcpp::Named("d") = d);
                        }
                        else
                        {
                            MoMALogger::error("Not implemented.");
                        }

                        five_d_list.insert(wrap_up, i, j, k, m, pc);

                        // Deflate the matrix
                        if (pc < rank - 1)
                        {
                            deflate();
                        }
                    }
                }
            }
        }
    }

    return five_d_list.get_list();
}

// 1. Return a list
// Rcpp::Named("lambda_u") = lambda_u,
// Rcpp::Named("lambda_v") = lambda_v,
// Rcpp::Named("alpha_u") = alpha_u,
// Rcpp::Named("alpha_v") = alpha_v,
// Rcpp::Named("u") = U,
// Rcpp::Named("v") = V,
// Rcpp::Named("d") = d
// 2. Dependence on MoMA's internal states: MoMA::X, MoMA::alpha_u/v, MoMA::lambda_u/v.
// 3. After calling MoMA::multi_rank, MoMA: MoMA::X becomes the corresponding deflated matrix.
// MoMA::u and MoMA::v become the leading penalized SVs of MoMA::X, using leading SVs of MoMA::X as
// start points.
Rcpp::List MoMA::multi_rank(int rank, arma::vec initial_u, arma::vec initial_v)
{
    if (rank <= 0)
    {
        MoMALogger::error("MoMA::multirank received non-positive rank: k = ") << rank;
    }
    // store results
    arma::mat U(X.n_rows, rank);
    arma::mat V(X.n_cols, rank);
    arma::vec d(rank);

    u = initial_u;
    v = initial_v;

    // find rank PCs
    for (int i = 0; i < rank; i++)
    {
        // Use MoMA::u and MoMA::v as start points.
        solve();
        U.col(i) = u;
        V.col(i) = v;
        d(i)     = arma::as_scalar(u.t() * X * v);
        // deflate X
        if (i < rank - 1)
        {
            deflate();
            // After deflation MoMA::u and MoMA::v are
            // re-initialized as leading SVs of the deflated matrix
            // MoMA::X = MoMA::X - d u v^T
        }
    }
    return Rcpp::List::create(Rcpp::Named("lambda_u") = lambda_u,
                              Rcpp::Named("lambda_v") = lambda_v, Rcpp::Named("alpha_u") = alpha_u,
                              Rcpp::Named("alpha_v") = alpha_v, Rcpp::Named("u") = U,
                              Rcpp::Named("v") = V, Rcpp::Named("d") = d);
}

// 1. Return a list
// Rcpp::Named("lambda_u") = lambda_u,
// Rcpp::Named("lambda_v") = lambda_v,
// Rcpp::Named("alpha_u") = alpha_u,
// Rcpp::Named("alpha_v") = alpha_v,
// Rcpp::Named("u") = U,
// Rcpp::Named("v") = V,
// Rcpp::Named("d") = d
// 2. Dependence on MoMA's internal states: MoMA::X.
// 3. After calling grid_search, MoMA::u and MoMA::v
// are the solution evaluated at the last grid point, using
// start points resulted by warm-start. MoMA::alpha_u/v, MoMA::lambda_u/v
// become the last grid point.
Rcpp::List MoMA::grid_search(const arma::vec &alpha_u,
                             const arma::vec &lambda_u,
                             const arma::vec &alpha_v,
                             const arma::vec &lambda_v,
                             arma::vec initial_u,
                             arma::vec initial_v)
{
    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u  = alpha_u.n_elem;
    int n_alpha_v  = alpha_v.n_elem;
    int n_total    = n_lambda_v * n_lambda_u * n_alpha_u * n_alpha_v;

    arma::mat U(X.n_rows, n_total);
    arma::mat V(X.n_cols, n_total);
    arma::vec d(n_total);

    int problem_id = 0;
    // MoMA::u and MoMA::v are initialzed as
    // the first leading SVs of MoMA::X
    for (int i = 0; i < n_lambda_u; i++)
    {
        for (int j = 0; j < n_lambda_v; j++)
        {
            for (int k = 0; k < n_alpha_u; k++)
            {
                for (int m = 0; m < n_alpha_v; m++)
                {
                    MoMALogger::info("Setting up model:")
                        << " lambda_u " << lambda_u(i) << " lambda_v " << lambda_v(j) << " alpha_u "
                        << alpha_u(k) << " alpha_v " << alpha_v(m);

                    set_penalty(lambda_u(i), lambda_v(j), alpha_u(k), alpha_v(m));

                    // MoMA::solve use the result from last
                    // iteration as starting point
                    solve();
                    U.col(problem_id) = u;
                    V.col(problem_id) = v;
                    d(problem_id)     = arma::as_scalar(u.t() * X * v);

                    problem_id++;
                }
            }
        }
    }
    if (problem_id != n_total)
    {
        MoMALogger::error("Internal error: solution not found for all grid points.");
    }
    return Rcpp::List::create(Rcpp::Named("lambda_u") = lambda_u,
                              Rcpp::Named("lambda_v") = lambda_v, Rcpp::Named("alpha_u") = alpha_u,
                              Rcpp::Named("alpha_v") = alpha_v, Rcpp::Named("u") = U,
                              Rcpp::Named("v") = V, Rcpp::Named("d") = d);
}
