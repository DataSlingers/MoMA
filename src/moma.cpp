// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "moma.h"
#include "moma_prox.h"
#include <algorithm>
#include <iostream>
#include <stdio.h>

enum class Solver{
    ISTA,
    FISTA
};


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double mat_norm(const arma::vec &u, const arma::mat &S_u)   // TODO: special case when S_u = I, i.e., alpha_u = 0.
{
    return arma::as_scalar(arma::sqrt(u.t() * S_u * u));
}


/////////////////
// Section 2: MoMA class
/////////////////

class MoMA{

private:
    /* matrix size */
    int n;
    int p;

    double prox_u_step_size; 
    double prox_v_step_size;
    double grad_u_step_size;
    double grad_v_step_size;
    
    Solver solver_type;
    arma::uword MAX_ITER;
    double EPS;

    const arma::mat &X; 
    //  careful about reference, if it refenrences 
    //  something that will be released in the 
    //  constructor, things go wrong
    
    // final results
    arma::vec u; 
    arma::vec v;
    // sparse penalty
    Prox *prox_u; // careful about memory leak and destructor stuff, can be replaced by Prox &prox_u;
    Prox *prox_v;
    // S = I + alpha*Omeg
    arma::mat S_u;  // to be special-cased
    arma::mat S_v;
    
    
public:
    ~MoMA(){
        delete prox_u;
        delete prox_v;
    }
    void check_valid();
    Solver string_to_SolverT(const std::string &s); // String to solver type {ISTA,FISTA}
    Prox* string_to_Proxptr(const std::string &s, double gamma);
   
    // turn user input into what we need to run the algorithm
    MoMA(const arma::mat &X_,   // note it is a reference
        /* sparsity*/
        std::string P_v,std::string P_u, 
        double lambda_v,double lambda_u,
        double gamma,
        /* smoothness */
        arma::mat Omega_u,arma::mat Omega_v,
        double alpha_u,double alpha_v,
        /* training para. */
        double i_EPS,arma::uword i_MAX_ITER,std::string i_solver):X(X_) // X has to be written in the initialization list
    {
        check_valid();
        MoMALogger::info("Setting up Our model\n");

       
        n = X.n_rows;
        p = X.n_cols;

        // Step 0: training para. setup
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
        grad_u_step_size = 1 / Lu;
        grad_v_step_size = 1 / Lv;
        prox_u_step_size = lambda_u / Lu;
        prox_v_step_size = lambda_v / Lv;

        // Step 2: initialize with SVD
        v = V.col(0);
        u = U.col(0);

        // Step 3: match proximal operator
        prox_u = string_to_Proxptr(P_u,gamma);
        prox_v = string_to_Proxptr(P_v,gamma);

        // Step 4: match gradient operator

    };

    void fit();
    Rcpp::List wrap(){ 
        u = u / norm(u);
        v = v / norm(v);
        double d = as_scalar(u.t() * X * v);
            return Rcpp::List::create(
            Rcpp::Named("u") = u,
            Rcpp::Named("v") = v,
            Rcpp::Named("d") = d,
            Rcpp::Named("DeflatedX") = X - d * u * v.t());
    }
};



void MoMA::check_valid(){
    MoMALogger::info("Checking input validity\n");
}

Solver MoMA::string_to_SolverT(const std::string &s){
    Solver res = Solver::ISTA;
    // we can first make s to upper case and provide more flexibility
    if (s.compare("ISTA") == 0)
        res = Solver::ISTA;
    else if (s.compare("FISTA") == 0)
        res = Solver::FISTA;
    else{
        MoMALogger::error("Your choice of algorithm not provided") << s;
    }
    return res;
}

Prox* MoMA::string_to_Proxptr(const std::string &s,double gamma){   // free it!
    Prox* res = nullptr;
    if (s.compare("LASSO") == 0)
       res = new Lasso();
    else if (s.compare("NONNEGLASSO") == 0)
        MoMALogger::error("Nonnegative Lasso not implemented yet!\n");
    else if (s.compare("SCAD") == 0)
        res = new SCAD(gamma);
    else if (s.compare("MCP") == 0)
        res = new MCP(gamma);
    else
        MoMALogger::error("Your sparse penalty is not provided!\n");
    return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List sfpca(
    const arma::mat &X ,

    arma::mat Omega_u,arma::mat Omega_v,  /* any idea to set up default values for these matrices? */
    double alpha_u = 0,double alpha_v = 0,

    std::string P_u = "LASSO",std::string P_v = "LASSO",
    double lambda_u = 0,double lambda_v = 0,
    double gamma = 3.7,

    double EPS = 1e-6,  
    long MAX_ITER = 1e+3,
    std::string solver = "ISTA"
)
{
    MoMA model(X,  
        /* sparsity*/
         P_v,P_u, 
        lambda_v,lambda_u,
        gamma,
        /* smoothness */
        Omega_u,Omega_v,
        alpha_u,alpha_v,
        /* optimizer parameter */
        EPS,
        MAX_ITER,
        solver);
    model.fit();

    return model.wrap();
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double test_norm(arma::vec x){
    return norm(x);
}

void MoMA::fit(){

        MoMALogger::info("Model info=========\n")<<"n:"<<n<<"\n"
            <<"p:" << p << "\n";
        MoMALogger::info("Start fitting.") << "\n";

        // keep the value of u at the start of outer loop, hence call it oldu1
        arma::vec oldu1 = arma::zeros<arma::vec>(n);
        arma::vec oldv1 = arma::zeros<arma::vec>(p);
        // keep the value of u at the start of inner loop
        arma::vec oldu2 = arma::zeros<arma::vec>(n);
        arma::vec oldv2 = arma::zeros<arma::vec>(p);

        // stopping tolerance
        int iter = 0;
        int iter_u = 0;
        int iter_v = 0;

        double in_u_tol = 1;   // tolerance for inner loop of u updates
        double in_v_tol = 1;   // tolerance for inner loop of v updates
        double out_tol = 1;    // that of outer loop

        if (solver_type == Solver::ISTA)
        {
            MoMALogger::info("Running ISTA!\n");
            MoMALogger::debug("==Before the loop: training setup==\n") 
                    << "\titer" << iter
                    << "\tEPS:" << EPS 
                    << "\tMAX_ITER:" << MAX_ITER << '\n';
            while (out_tol > EPS && iter < MAX_ITER)
            {
                
                oldu1 = u;  
                oldv1 = v;
                in_u_tol = 1;
                in_v_tol = 1;
                iter_u = 0;
                iter_v = 0;
                while (in_u_tol > EPS)
                {
                    iter_u++;
                    oldu2 = u;  
                    // gradient step
                    u = u + grad_u_step_size * (X*v - S_u*u);  // TODO: special case when alpha_u = 0 => S_u = I
                    // proximal step
                    u = prox_u->prox(u,prox_u_step_size);
                    // nomalize w.r.t S_u
                    norm(u) > 0 ? u /= mat_norm(u, S_u) : u.zeros();    // Sometimes mat_norm(u,S_u) is so close to zero that u becomes NaN

                    in_u_tol = norm(u - oldu2) / norm(oldu2);
                //    if(iter_u %100 == 0)
                        MoMALogger::debug("---update u ") << iter_u << "--\n" 
                            << "in_u_tol:" << in_u_tol << "\t iter" << iter_u;
                }

                while (in_v_tol > EPS)
                {
                    iter_v++;
                    oldv2 = v;
                    // gradient step
                    v = v + grad_v_step_size * (X.t()*u - S_v*v);    // TODO: special case
                    // proximal step
                    v = prox_v->prox(v,prox_v_step_size);
                    norm(v) > 0 ? v /= mat_norm(v, S_v) : v.zeros();
                    in_v_tol = norm(v - oldv2) / norm(oldv2);
                    // if(iter_v %100 == 0)
                        MoMALogger::debug("---update v ") << iter_v << "---\n"
                            << "in_v_tol:" << in_v_tol << "\t iter" << iter_v;
                }

                out_tol = norm(oldu1 - u) / norm(oldu1) + norm(oldv1 - v) / norm(oldv1);
                iter++;
                MoMALogger::debug("--Finish iter:") << iter << "---\n";
            }
        }
        else if (solver_type == Solver::FISTA){
            MoMALogger::error("FISTA is not provided yet!\n");
        }
        else{
            MoMALogger::error("Your choice of solver is not provided yet!");
        }
        MoMALogger::debug("==After the outer loop!==\n") 
                   << "out_tol:" << out_tol << "\t iter" << iter;
    }