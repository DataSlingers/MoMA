# include "moma_prox.h"

// [[Rcpp::export]]
arma::vec test_prox_lasso(const arma::vec &x, double l)
{
    Lasso a;
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_nnlasso(const arma::vec &x, double l)
{
    NonNegativeLasso a;
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_scad(const arma::vec &x, double l, double gamma = 3.7)
{
    SCAD a(gamma);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_scadvec(const arma::vec &x, double l, double gamma = 3.7)
{

    SCAD a(gamma);
    return a.vec_prox(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_nnscad(const arma::vec &x, double l, double gamma = 3.7)
{
    NonNegativeSCAD a(gamma);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_mcp(const arma::vec &x, double l, double gamma = 4)
{
    MCP a(gamma);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_mcpvec(const arma::vec &x, double l, double gamma = 4)
{
    MCP a(gamma);
    return a.vec_prox(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_nnmcp(const arma::vec &x, double l, double gamma = 4)
{
    NonNegativeMCP a(gamma);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_grplasso(const arma::vec &x, const arma::vec &g, double l)
{
    GrpLasso a(g);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_nngrplasso(const arma::vec &x, const arma::vec &g, double l)
{
    NonNegativeGrpLasso a(g);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_orderedfusion(const arma::vec &x,double l)
{
    OrderedFusedLasso a;
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_spfusedlasso(const arma::vec &x,double l,double lambda2)
{
    // lambda2: the level of penalty on
    // the absolute values of the coefficients
    SparseFusedLasso a(lambda2);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_fusion(const arma::vec &x,double l,const arma::mat w,bool ADMM,bool acc,double prox_eps=1e-10)
{
    Fusion a(w,ADMM,acc,prox_eps);
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_l1gf(const arma::vec &x,double l,int k = 1)
{
    L1TrendFiltering a(x.n_elem,k);
    return a(x,l);
}

// [[Rcpp::export]]
int test_df_orderedfusion(const arma::vec &x)
{
    OrderedFusedLasso a;
    return a.df(x);
}

// [[Rcpp::export]]
int test_df_spfusedlasso(const arma::vec &x,double lambda2)
{
    // lambda2: the level of penalty on
    // the absolute values of the coefficients
    SparseFusedLasso a(lambda2);
    return a.df(x);
}

// [[Rcpp::export]]
int test_df_l1gf(const arma::vec &x,int k = 1)
{
    L1TrendFiltering a(x.n_elem,k);
    return a.df(x);
}

// [[Rcpp::export]]
int test_df_grplasso(const arma::vec &x, const arma::vec &g)
{
    GrpLasso a(g);
    return a.df(x);
}
