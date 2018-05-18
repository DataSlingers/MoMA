// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "moma.h"
#include "moma_prox.h"
#include <algorithm>
using namespace Rcpp;
using namespace arma;
using namespace std;
enum class Solver{
    ISTA,
    FISTA
};
enum class SparsityType{
    LASSO,
    NONEGLASSO,
    SCAD,
    MCP
};

/////////////////
// Section 2: MoMA class
/////////////////



class MoMA{

private:
    // Input parameters
    double prox_u_step; //lambda_u, L_u(S_u(Omeg,alpha_a)
    double prox_v_step;
    double grad_u_step;
    double grad_v_step;
    
    Solver solver_type;

    long MAX_ITER;
    double EPS;

    // final results
    arma::vec u; 
    arma::vec v;
    // Model, sparse penalty, non_neg
    Prox *prox_u; // careful about release space and destructor stuff
    Prox *prox_v;
    // alpha=0, S_u = I
    arma::mat S_u;  // to be special case
    arma::mat S_v;
    
public:
    // turn user input into what we need to run the algorithm
    MoMA(arma::mat X,
        // sparsity
        std::string P_v,
        std::string P_u, // assume for now they have same type of penalty
        double lambda_v,
        double lambda_u,
        double gamma,
        bool non_neg, // 1 means activate non-negativity constraint
        // smoothness
        arma::mat Omega_u,
        arma::mat Omega_v,
        double alpha_u,
        double alpha_v,
        // training para.
        double i_EPS,
        long i_MAX_ITER,
        std::string i_solver)
    {
    check_valid();
    Rcpp::Rcout<< "Setting up\n";

    int n = X.n_rows, p = X.n_cols;
    // Step 0: easy setup
    MAX_ITER = i_MAX_ITER;
    EPS = i_EPS;
    solver_type = string_to_ST(i_solver);

    // Step 1: find Su,Sv, and thus Lu,Lv (largest eigenvalue), and thus stepsizes
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    S_u.eye(arma::size(Omega_u));
    S_v.eye(arma::size(Omega_v));
    S_u += n * alpha_u * Omega_u;
    S_v += p * alpha_v * Omega_v;

    double Lu = arma::eig_sym(S_u).max() + 0.01; // +0.01 for convergence
    double Lv = arma::eig_sym(S_v).max() + 0.01;
    
    grad_u_step = 1 / Lu;
    grad_v_step = 1 / Lv;
    prox_u_step = lambda_u / Lu;
    prox_v_step = lambda_v / Lv;

    // Step 2: initialize with SVD
    v = V.col(0);
    u = U.col(0);

    };
    void check_valid();


    Solver string_to_ST(const std::string &s){
        // we can first make s to upper case and provide more flexibility
        if (s.compare("ISTA") == 0)
            return Solver::ISTA;
        else if (s.compare("FISTA") == 0)
            return Solver::FISTA;
        else
            MoMALogger::error("Your choice of algorithm not provided") << s;
    }
};



void MoMA::check_valid(){
    Rcpp::Rcout << "Checking input validity\n";
}