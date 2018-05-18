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

// // 
// class Grad{
// protected:
//     arma::vec b;  // keeps changing when training
// public:
//     Grad() = delete; // has to be initialized
//     Grad(const arma::vec &b_){b=b_;}
//     void update(const arma::vec &b_){b=b_;}
//     virtual arma::vec desc(const arma::vec &x, double stepsize){
//         return x - stepsize * (x + b);
//     }
// };

// class Comm:Grad{    // Common case when A != I
//     arma::mat A;    // unchanged during training
// public:
//     Comm() = delete;
//     Comm(const arma::vec &b_,const arma::mat &A_){
//         A = A_;
//         b = b_;
//     }
//     arma::vec desc(const arma::vec &x, double stepsize){
//         return x - stepsize * (A*x + b);
//     }
// };


class MoMA{

private:
    
    double prox_u_step; 
    double prox_v_step;
    double grad_u_step;
    double grad_v_step;
    
    Solver solver_type;
    long MAX_ITER;
    double EPS;

    // final results
    arma::vec u; 
    arma::vec v;
    // sparse penalty
    Prox *prox_u; // careful about release space and destructor stuff
    Prox *prox_v;
    // S = I + alpha*Omeg
    arma::mat S_u;  // to be special case
    arma::mat S_v;
    
    
public:

    void check_valid();
    Solver string_to_SolverT(const std::string &s);
    Prox* string_to_SparseT(const std::string &s,double gamma);
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
        std::string i_solver){
    check_valid();
    Rcpp::Rcout<< "Setting up\n";

    int n = X.n_rows, p = X.n_cols;
    // Step 0: easy setup
    MAX_ITER = i_MAX_ITER;
    EPS = i_EPS;
    solver_type = string_to_SolverT(i_solver);

    // Step 1: find Su,Sv
    arma::mat U;
    arma::vec s;
    arma::mat V;
    arma::svd(U, s, V, X);
    S_u.eye(arma::size(Omega_u));
    S_v.eye(arma::size(Omega_v));
    S_u += n * alpha_u * Omega_u;
    S_v += p * alpha_v * Omega_v;

    // Step 1.2: find Lu,Lv
    double Lu = arma::eig_sym(S_u).max() + 0.01; // +0.01 for convergence
    double Lv = arma::eig_sym(S_v).max() + 0.01;
    
    // Step 1.3: all kinds of stepsize
    grad_u_step = 1 / Lu;
    grad_v_step = 1 / Lv;
    prox_u_step = lambda_u / Lu;
    prox_v_step = lambda_v / Lv;

    // Step 2: initialize with SVD
    v = V.col(0);
    u = U.col(0);

    // Step 3: match proximal operator
    prox_u = string_to_SparseT(P_u,gamma);
    prox_v = string_to_SparseT(P_v,gamma);

    // Step 4: match gradient operator

    };



};



void MoMA::check_valid(){
    Rcpp::Rcout << "Checking input validity\n";
}

Solver MoMA::string_to_SolverT(const std::string &s){
    // we can first make s to upper case and provide more flexibility
    if (s.compare("ISTA") == 0)
        return Solver::ISTA;
    else if (s.compare("FISTA") == 0)
        return Solver::FISTA;
    // // assume `s` is in the set {ISTA, FISTA}
    // else
    //     MoMALogger::error("Your choice of algorithm not provided") << s;
}

Prox* MoMA::string_to_SparseT(const std::string &s,double gamma){
    if (s.compare("LASSO") == 0)
       return new Lasso();
    else if (s.compare("NONNEGLASSO") == 0)
        MoMALogger::error("Nonnegative Lasso not implemented yet!\n");
    else if (s.compare("SCAD") == 0)
        return new Scad(gamma);
    else if (s.compare("MCP") == 0)
        return new Mcp(gamma);
    // // assume `s` is in the set {ISTA, FISTA}
    // else
    //     MoMALogger::error("Your sparse penalty is not provided!\n");
}