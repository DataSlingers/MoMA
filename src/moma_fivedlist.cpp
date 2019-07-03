#include "moma_fivedlist.h"

RcppFiveDList::RcppFiveDList(int n_alpha_u, int n_lambda_u, int n_alpha_v, int n_lambda_v, int k)
    : n_alpha_u(n_alpha_u),
      n_lambda_u(n_lambda_u),
      n_alpha_v(n_alpha_v),
      n_lambda_v(n_lambda_v),
      k(k),
      flattened_list(n_alpha_u * n_alpha_v * n_lambda_u * n_lambda_v * k)
{
    flattened_list.attr("dim") =
        Rcpp::NumericVector::create(n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v, k);
    flattened_list.attr("class") = "MoMA_5D_list";
};

int RcppFiveDList::insert(Rcpp::List object,
                          int alpha_u_i,
                          int lambda_u_i,
                          int alpha_v_i,
                          int lambda_v_i,
                          int k_i)
{
    // insert object in the alpha_u_i-th position along the alpha_u-axis
    // and so on
    if (alpha_u_i < 0 || alpha_u_i >= n_alpha_u || lambda_u_i < 0 || lambda_u_i >= n_lambda_u ||
        alpha_v_i < 0 || alpha_v_i >= n_alpha_v || lambda_v_i < 0 || lambda_v_i >= n_lambda_v ||
        k_i < 0 || k_i >= k)
    {
        MoMALogger::error("Invalid index is passed to RcppFiveDList::insert. ")
            << "Dimension is (" << n_alpha_u << ", " << n_lambda_u << ", " << n_alpha_v << ", "
            << n_lambda_v << "," << k << "), received (" << alpha_u_i << ", " << lambda_u_i << ", "
            << alpha_v_i << ", " << lambda_v_i << ", " << k_i << ").";
    }

    flattened_list(k_i + k * (lambda_v_i +
                              n_lambda_v * (alpha_v_i +
                                            n_alpha_v * (lambda_u_i + n_lambda_u * (alpha_u_i))))) =
        object;
    return 0;
}

Rcpp::List RcppFiveDList::get_list()
{
    return flattened_list;
}
