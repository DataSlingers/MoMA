# include "moma_prox.h"
# include "moma_solver.h"
# include "moma.h"

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
arma::vec test_prox_fusedlassopath(const arma::vec &x,double l)
{
    OrderedFusedLasso a;
    return a(x,l);
}

// [[Rcpp::export]]
arma::vec test_prox_fusedlassodp(const arma::vec &x,double l)
{
    OrderedFusedLassoDP a;
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
arma::vec test_prox_slope(const arma::vec &x,double l)
{
    // lambda2: the level of penalty on
    // the absolute values of the coefficients
    SLOPE a(x.n_elem);
    return a(x,l);
}

// [[Rcpp::export]]
int test_df_orderedfusion(const arma::vec &x)
{
    OrderedFusedLasso a;
    return a.df(x);
}

// [[Rcpp::export]]
int test_df_spfusedlasso(const arma::vec &x)
{
    // SparseFusedLasso object needs a `lambda2` argument
    // (the level of penalty on the absolute values of
    // the coefficients) to initialize
    SparseFusedLasso a(1);
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

// [[Rcpp::export]]
double test_BIC( 
        const arma::vec y, const arma::vec y_est,
        const std::string &algorithm_string,
        double i_alpha, const arma::mat &i_Omega,
        double i_lambda, Rcpp::List prox_arg_list,
        int dim,
        double i_EPS=1e-6, int i_MAX_ITER=1e+3)
{
        
    PR_solver solver(
            algorithm_string,
            i_alpha, i_Omega,
            i_lambda, prox_arg_list,
            i_EPS, i_MAX_ITER, dim);

    return solver.bic(y, y_est);
}

// This function solves a squence of lambda's and alpha's
// [[Rcpp::export]]
Rcpp::List testnestedBIC(
    const arma::mat &X,    // We should not change any variable in R, so const ref
    const arma::vec &alpha_u,
    const arma::vec &alpha_v,
    const arma::mat &Omega_u, // Default values for these matrices should be set in R
    const arma::mat &Omega_v,
    const arma::vec &lambda_u,
    const arma::vec &lambda_v,
    const Rcpp::List &prox_arg_list_u,
    const Rcpp::List &prox_arg_list_v,
    double EPS,
    long MAX_ITER,
    double EPS_inner,
    long MAX_ITER_inner,
    std::string solver,
    int selection_criterion_alpha_u,  // flags; = 0 means grid, = 1 means BIC search
    int selection_criterion_alpha_v,
    int selection_criterion_lambda_u,
    int selection_criterion_lambda_v,
    int k = 1){


    int n_lambda_u = lambda_u.n_elem;
    int n_lambda_v = lambda_v.n_elem;
    int n_alpha_u = alpha_u.n_elem;
    int n_alpha_v = alpha_v.n_elem;

    if(n_lambda_v == 0 || n_lambda_u == 0 || n_alpha_u == 0 || n_alpha_v == 0){
        MoMALogger::error("Please specify all four parameters.");
    }

    int n_total = n_lambda_v * n_lambda_u * n_alpha_u * n_alpha_v;

    // NOTE: arguments should be listed
    // in the exact order of MoMA constructor
    MoMA problem(X,
              /* sparsity */
              lambda_u(0),
              lambda_v(0),
              prox_arg_list_u,
              prox_arg_list_v,
              /* smoothness */
              alpha_u(0),
              alpha_v(0),
              Omega_u,
              Omega_v,
              /* algorithm parameters */
              EPS,
              MAX_ITER,
              EPS_inner,
              MAX_ITER_inner,
              solver);
    return problem.grid_BIC_mix(alpha_u,alpha_v,lambda_u,lambda_v,
                        selection_criterion_alpha_u,selection_criterion_alpha_v,selection_criterion_lambda_u,selection_criterion_lambda_v);
    // return problem.select_nestedBIC(alpha_u,alpha_v,lambda_u,lambda_v,5);
}
