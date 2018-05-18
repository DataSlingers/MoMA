// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma.h"
#include "moma_prox.h"
using namespace Rcpp;
using namespace arma;
using namespace std;
typedef enum { ISTA, FISTA } SOLVER_TYPE;


/////////////////
// Section 2: MoMA class
/////////////////
std::string str_toupper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), 
                   [](unsigned char c){ return std::toupper(c); });
    return s;
}

class Grad{  // In our case, we only have grad of such form: f'(x) = Ax + b
public:
    arma::mat A;    // remains unchanged
    arma::vec b;    // keeps changing during training
};

class MoMA{

private:
    // Input parameters
  
    double prox_u_step; //lambda_u, L_u(S_u(Omeg,alpha_a)
    double prox_v_step;
    double grad_u_step;
    double grad_v_step;
    
    SOLVER_TYPE solver_type;

    long MAX_ITER;
    double EPS;

    // final results
    arma::vec u; 
    arma::vec v;
    // Model, sparse penalty, non_neg
    Prox *prox_u; 
    Prox *prox_v;
    // alpha ==0, S_u = I
    Grad *grad_u; 
    Grad *grad_v;
    arma::mat Su;
    arma::mat Sv;

public:
    // turn user input into what we need to run the algorithm
    void setup(arma::mat X,
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
        double EPS,
        long MAX_ITER,
        std::string solver_type_string
    ){
        check_valid();
        Rcpp::Rcout<< "Converting\n";
    };
    void check_valid();
};



void MoMA::check_valid(){
    Rcpp::Rcout << "Checking input validity\n";
}