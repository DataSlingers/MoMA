#ifndef moma_fivedlist_H
#define moma_fivedlist_H 1
#include "moma.h"

class RcppFiveDList
{
    int n_alpha_u;
    int n_lambda_u;
    int n_alpha_v;
    int n_lambda_v;
    int k;
    Rcpp::List flattened_list;

  public:
    RcppFiveDList(int n_alpha_u, int n_lambda_u, int n_alpha_v, int n_lambda_v, int k = 1);

    int insert(Rcpp::List object,
               int alpha_u_i,
               int lambda_u_i,
               int alpha_v_i,
               int lambda_v_i,
               int k_i = 0);

    Rcpp::List get_list();
};

#endif
