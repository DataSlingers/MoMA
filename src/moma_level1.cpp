// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil;
// -*-
#include "moma.h"

// auxiliary functions for MoMA::grid_BIC_mix
const arma::vec &set_greedy_grid(const arma::vec &grid, int want_grid)
{
    if (want_grid == 1)
    {
        return grid;
    }
    else if (want_grid == 0)
    {
        return MOMA_EMPTY_GRID_OF_LENGTH1;
    }
}

arma::vec set_bic_grid(const arma::vec &grid, int want_bic, int i)
{
    if (want_bic == 1)
    {
        return grid;
    }
    else if (want_bic == 0)
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

        double scale_u = norm(oldu) == 0.0 ? 1 : norm(oldu);
        double scale_v = norm(oldv) == 0.0 ? 1 : norm(oldv);

        tol = norm(oldu - u) / scale_u + norm(oldv - v) / scale_v;
        MoMALogger::debug("Finish nested greedy BIC search outer loop. (iter, tol) = (")
            << iter << "," << tol << "), "
            << "(bic_u, bic_v) = (" << (double)u_result["bic"] << "," << (double)v_result["bic"]
            << ")";
    }

    // A final run on the selected parameter
    if (final_run)
    {
        MoMALogger::message("Start a final run on the chosen parameters.");
        double opt_lambda_u = u_result["lambda"];
        double opt_lambda_v = v_result["lambda"];
        double opt_alpha_u  = u_result["alpha"];
        double opt_alpha_v  = v_result["alpha"];

        reset(opt_lambda_u, opt_lambda_v, opt_alpha_u, opt_alpha_v);
        initialize_uv();
        solve();  // Use MoMA::u and MoMA::v as start points

        u_result["vector"] = u;
        v_result["vector"] = v;
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
                              int selection_criterion_alpha_u,  // flags; = 0 means grid, =
                                                                // 01 means BIC search
                              int selection_criterion_alpha_v,
                              int selection_criterion_lambda_u,
                              int selection_criterion_lambda_v,
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
    const arma::vec &grid_lu = set_greedy_grid(lambda_u, !selection_criterion_lambda_u);
    const arma::vec &grid_lv = set_greedy_grid(lambda_v, !selection_criterion_lambda_v);
    const arma::vec &grid_au = set_greedy_grid(alpha_u, !selection_criterion_alpha_u);
    const arma::vec &grid_av = set_greedy_grid(alpha_v, !selection_criterion_alpha_v);

    // Test that if a grid is set to be BIC-search grid, then
    // the above code should set grid_xx to the vector [-1]
    if ((selection_criterion_alpha_u == 1 && (grid_au.n_elem != 1 || grid_au(0) != -1)) ||
        (selection_criterion_alpha_v == 1 && (grid_av.n_elem != 1 || grid_av(0) != -1)) ||
        (selection_criterion_lambda_u == 1 && (grid_lu.n_elem != 1 || grid_lu(0) != -1)) ||
        (selection_criterion_lambda_v == 1 && (grid_lv.n_elem != 1 || grid_lv(0) != -1)))
    {
        MoMALogger::error("Wrong grid-search grid!")
            << "grid_lu.n_elem=" << grid_lu.n_elem << ", grid_av.n_elem=" << grid_av.n_elem
            << ", grid_lu.n_elem" << grid_lu.n_elem << ", grid_lv.n_elem" << grid_lv.n_elem;
    }

    int n_lambda_u = grid_lu.n_elem;
    int n_lambda_v = grid_lv.n_elem;
    int n_alpha_u  = grid_au.n_elem;
    int n_alpha_v  = grid_av.n_elem;

    RcppFiveDList five_d_list(n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v, rank);

    arma::mat original_X =
        X;  // keep a copy of X, because MoMA::X will be contanminated during defaltion
    for (int i = 0; i < n_alpha_u; i++)
    {
        for (int j = 0; j < n_lambda_u; j++)
        {
            for (int k = 0; k < n_alpha_v; k++)
            {
                for (int m = 0; m < n_lambda_v; m++)
                {
                    arma::vec bic_au_grid = set_bic_grid(alpha_u, selection_criterion_alpha_u, i);
                    arma::vec bic_lu_grid = set_bic_grid(lambda_u, selection_criterion_lambda_u, j);
                    arma::vec bic_av_grid = set_bic_grid(alpha_v, selection_criterion_alpha_v, k);
                    arma::vec bic_lv_grid = set_bic_grid(lambda_v, selection_criterion_lambda_v, m);

                    if ((selection_criterion_alpha_u == 0 && bic_au_grid.n_elem != 1) ||
                        (selection_criterion_alpha_v == 0 && bic_av_grid.n_elem != 1) ||
                        (selection_criterion_lambda_u == 0 && bic_lu_grid.n_elem != 1) ||
                        (selection_criterion_lambda_v == 0 && bic_lv_grid.n_elem != 1))
                    {
                        MoMALogger::error("Wrong BIC search grid!");
                    }

                    set_X(original_X);
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

                        five_d_list.insert(
                            Rcpp::List::create(Rcpp::Named("u") = u_result,
                                               Rcpp::Named("v") = v_result, Rcpp::Named("k") = pc,
                                               Rcpp::Named("X") = X),
                            i, j, k, m, pc);

                        // Deflate the matrix
                        if (pc < rank - 1)
                        {
                            arma::vec curu = Rcpp::as<Rcpp::NumericVector>(u_result["vector"]);
                            arma::vec curv = Rcpp::as<Rcpp::NumericVector>(v_result["vector"]);
                            double d       = arma::as_scalar(curu.t() * X * curv);
                            deflate(d);
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
            deflate(d(i));
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

                    reset(lambda_u(i), lambda_v(j), alpha_u(k), alpha_v(m));

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