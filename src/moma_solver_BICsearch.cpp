#include "moma_solver_BICsearch.h"

void BIC_searcher::bind(PR_solver *object, Criterion method)
{
    pr_solver = object;
    cri       = method;
}

double BIC_searcher::cur_criterion(arma::vec y, const arma::vec &est)
{
    return (pr_solver->*cri)(y, est);
}

// Return a Rcpp::List:
//   Rcpp::Named("lambda") = opt_lambda_u,
//   Rcpp::Named("alpha")  = opt_alpha_u,
//   Rcpp::Named("vector") = working_selected_u,
//   Rcpp::Named("bic")    = minbic_u
Rcpp::List BIC_searcher::search(const arma::vec &y,          // min_{u} || y - u || + ...penalty...
                                const arma::vec &initial_u,  // start point
                                const arma::vec &alpha_u,
                                const arma::vec &lambda_u)
{
    arma::vec working_selected_u;
    arma::vec working_u = initial_u;
    double working_bic_u;
    double minbic_u = MOMA_INFTY;
    double opt_alpha_u;
    double opt_lambda_u;
    for (int i = 0; i < alpha_u.n_elem; i++)
    {
        // Put lambda_u in the inner loop to avoid reconstructing S many times
        for (int j = 0; j < lambda_u.n_elem; j++)
        {
            pr_solver->set_penalty(lambda_u(j), alpha_u(i));
            // working_u is the solution of the previous problem
            working_u     = pr_solver->solve(y, working_u);
            working_bic_u = cur_criterion(y, working_u);
            MoMALogger::debug("(curBIC, minBIC, lambda, alpha) = (")
                << working_bic_u << "," << minbic_u << "," << lambda_u(j) << "," << alpha_u(i)
                << ")";
            if (working_bic_u < minbic_u)
            {
                minbic_u           = working_bic_u;
                working_selected_u = working_u;
                opt_lambda_u       = lambda_u(j);
                opt_alpha_u        = alpha_u(i);
            }
        }
    }
    MoMALogger::debug("Finish greedy BIC, chosen (minBIC, alpha, lambda) = (")
        << minbic_u << ", " << opt_alpha_u << ", " << opt_lambda_u << ").";

    return Rcpp::List::create(
        Rcpp::Named("lambda") = opt_lambda_u, Rcpp::Named("alpha") = opt_alpha_u,
        Rcpp::Named("vector") = working_selected_u, Rcpp::Named("bic") = minbic_u);
};
