// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma.h"

enum class Solver{
    ISTA,
    FISTA
};

inline double mat_norm(const arma::vec &u, const arma::mat &S_u){
    // TODO: special case when S_u = I, i.e., alpha_u = 0.
    return arma::as_scalar(arma::sqrt(u.t() * S_u * u));
}

class MoMA{

private:
    /* matrix size */
    int n; // rows
    int p; // columns

    // Step size for proximal gradient algorithm
    //   - since this is a linear model internally, we can used a fixed
    //     step size without backtracking
    double prox_u_step_size;
    double prox_v_step_size;
    double grad_u_step_size;
    double grad_v_step_size;

    Solver solver_type;
    arma::uword MAX_ITER;
    double EPS;

    // Be careful when using references: if it references
    // something that is released during the call to the constructor,
    // things can go wrong
    const arma::mat& X;

    // Results -- will be modified during iterations and copied back to R
    arma::vec u;
    arma::vec v;

    // Proximal operators for sparsity inducing penalties
    //
    // Note that currently the threshold level is not defined in the Prox object
    //
    // We are using raw pointers here, so need to be careful to free them in ~MoMA()
    // TODO: Use const references for the prox objects
    Prox *prox_u;
    Prox *prox_v;

    // S = I + alpha * Omega for u, v smoothing
    // TODO: Special-case the Omega = 0 => S = I case
    arma::mat S_u;
    arma::mat S_v;


public:
    ~MoMA(){
        delete prox_u;
        delete prox_v;
    }

    void check_valid();
    Solver string_to_SolverT(const std::string &s); // String to solver type {ISTA,FISTA}
    Prox* string_to_Proxptr(const std::string &s,double gamma,const arma::vec &group,bool nonneg);
 
    // Parse user input into a MoMA object which defines the problem and algorithm
    // used to solve it.
    //
    // TODO: Decouple problem defintion and algorithmic choices
    //
    MoMA(const arma::mat &X_, // Pass X_ as a reference to avoid copy
        /*
         * sparsity - enforced through penalties
         */
        std::string P_v, // Sparsity penalty info
        std::string P_u,
        double lambda_v, // regularization level
        double lambda_u,
        double gamma,    // Non-convexity parameter
        bool nonneg_u,   // Non-negativity indicator
        bool nonneg_v,
        /* 
        * grouping 
        */
        const arma::vec &group_u,
        const arma::vec &group_v,
        /*
         * smoothness - enforced through constraints
         */
        arma::mat Omega_u, // Smoothing matrices
        arma::mat Omega_v,
        double alpha_u,    // Smoothing levels
        double alpha_v,
        /*
         * Algorithm parameters:
         */
        double i_EPS,
        arma::uword i_MAX_ITER,
        std::string i_solver): X(X_) // const reference must be passed to initializer list
    {
        check_valid();
        MoMALogger::info("Setting up our model");

        n = X.n_rows;
        p = X.n_cols;

        // Step 0: Set optimizer parameters
        MAX_ITER = i_MAX_ITER;
        EPS = i_EPS;
        solver_type = string_to_SolverT(i_solver);

        // Step 1a: Calculate smoothing matrices
        arma::mat U;
        arma::vec s;
        arma::mat V;
        arma::svd(U, s, V, X);
        S_u.eye(arma::size(Omega_u));
        S_v.eye(arma::size(Omega_v));
        S_u += n * alpha_u * Omega_u;
        S_v += p * alpha_v * Omega_v;

        // Step 1b: Calculate leading eigenvalues of smoothing matrices
        //          -> used for prox gradient step sizes
        double Lu = arma::eig_sym(S_u).max() + MOMA_EIGENVALUE_REGULARIZATION;
        double Lv = arma::eig_sym(S_v).max() + MOMA_EIGENVALUE_REGULARIZATION;

        // Step 1c: Set step sizes
        grad_u_step_size = 1 / Lu;
        grad_v_step_size = 1 / Lv;
        prox_u_step_size = lambda_u / Lu;
        prox_v_step_size = lambda_v / Lv;

        // Step 2: Initialize to leading singular vectors
        //
        //         MoMA is a regularized SVD, which is a non-convex (bi-convex)
        //         problem, so we need to be cautious about initialization to
        //         avoid local-minima. Initialization at the SVD (global solution
        //         to the non-regularized problem) seems to be a good trade-off:
        //         for problems with little regularization, the MoMA solution will
        //         lie near the SVD solution; for problems with significant regularization
        //         the problem becomes more well-behaved and less sensitive to
        //         initialization
        v = V.col(0);
        u = U.col(0);

        // Step 3: Construct proximal operators
        prox_u = string_to_Proxptr(P_u,gamma,group_u,nonneg_u);
        prox_v = string_to_Proxptr(P_v,gamma,group_v,nonneg_v);
    };

    void fit(); // Implemented below

    Rcpp::List wrap(){
        // Wrap results before returning to R

        u = u / norm(u); // Normalize one more time just in case
        v = v / norm(v);
        double d = as_scalar(u.t() * X * v); // Calculate singular value

        return Rcpp::List::create(
                     Rcpp::Named("u") = u,
                     Rcpp::Named("v") = v,
                     Rcpp::Named("d") = d,
                     Rcpp::Named("DeflatedX") = X - d * u * v.t());
    }
};

void MoMA::check_valid(){
    // TODO -- Check input
    MoMALogger::info("Checking input validity");
}

Solver MoMA::string_to_SolverT(const std::string &s){
    Solver res = Solver::ISTA;
    // TODO: capitalize s to be more robust to user-error
    if (s.compare("ISTA") == 0)
        res = Solver::ISTA;
    else if (s.compare("FISTA") == 0)
        res = Solver::FISTA;
    else{
        MoMALogger::error("Your choice of algorithm not provided") << s;
    }
    return res;
}

Prox* MoMA::string_to_Proxptr(const std::string &s,double gamma,const arma::vec &group,bool nonneg){
    // IMPORTANT: this must be freed somewhere
    Prox* res = new Prox();
    if (s.compare("LASSO") == 0){
        if(nonneg)
            res = new NonNegativeLasso();
        else
            res = new Lasso();
    }
    else if (s.compare("SCAD") == 0){
        if(nonneg)
            res = new NonNegativeSCAD(gamma);
        else 
            res = new SCAD(gamma);
    }
    else if (s.compare("MCP") == 0){
        if(nonneg)   
            res = new NonNegativeMCP(gamma);
        else 
            res = new MCP(gamma);
    }
    else if(s.compare("GRPLASSO") == 0){
        if(nonneg)
            res = new NonNegativeGrpLasso(group);
        else
            res = new GrpLasso(group);
    }
    else if(s.compare("FUSION") == 0){
           MoMALogger::error("Fusion is not provided yet");
    }
    else
        MoMALogger::warning("Your sparse penalty is not provided by us/specified by you! Use `Prox` by default");
    return res;
}

// [[Rcpp::export]]
Rcpp::List sfpca(
    const arma::mat& X,
    arma::mat Omega_u, // Default values for these matrices should be set in R
    arma::mat Omega_v,
    double alpha_u = 0,
    double alpha_v = 0,
    std::string P_u = "LASSO",
    std::string P_v = "LASSO",
    double lambda_u = 0,
    double lambda_v = 0,
    double gamma = 3.7,
    bool nonneg_u = 0, 
    bool nonneg_v = 0,
    arma::vec group_u = Rcpp::IntegerVector::create(0), 
    arma::vec group_v = Rcpp::IntegerVector::create(0),
    double EPS = 1e-6,
    long MAX_ITER = 1e+3,
    std::string solver = "ISTA"){

    MoMA model(X,
              /* sparsity */
              P_v,
              P_u,
              lambda_v,
              lambda_u,
              gamma,
              /* non-negativity */
              nonneg_u, 
              nonneg_v,
              /* grouping */
              group_u,
              group_v,
              /* smoothness */
              Omega_u,
              Omega_v,
              alpha_u,
              alpha_v,
              /* algorithm parameters */
              EPS,
              MAX_ITER,
              solver);

    model.fit(); // Run MoMA!

    return model.wrap();
}

void MoMA::fit(){
    MoMALogger::info("Model info=========\n") << "n:" << n <<"\n"
                                              << "p:" << p << "\n";
    MoMALogger::info("Start fitting.");

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

    if (solver_type == Solver::ISTA){
        MoMALogger::info("Running ISTA!");
        MoMALogger::debug("==Before the loop: training setup==\n")
                    << "\titer" << iter
                    << "\tEPS:" << EPS
                    << "\tMAX_ITER:" << MAX_ITER << '\n';
        while (out_tol > EPS && iter < MAX_ITER) {
            oldu1 = u;
            oldv1 = v;
            in_u_tol = 1;
            in_v_tol = 1;
            iter_u = 0;
            iter_v = 0;
            while (in_u_tol > EPS){
                iter_u++;
                oldu2 = u;
                // Gradient step
                // TODO: special case when alpha_u = 0 => S_u = I
                u = u + grad_u_step_size * (X*v - S_u*u);

                // Proximal step
                u = (*prox_u)(u,prox_u_step_size);
                // Normalize with respect to S_u
                // Sometimes mat_norm(u, S_u) is so close to zero that u becomes NaN
                norm(u) > 0 ? u /= mat_norm(u, S_u) : u.zeros();

                in_u_tol = norm(u - oldu2) / norm(oldu2);
                //    if(iter_u % 100 == 0)
                        MoMALogger::debug("---update u ") << iter_u << "--\n"
                                                          << "in_u_tol:" << in_u_tol
                                                          << "\t iter" << iter_u;
            }

            while (in_v_tol > EPS) {
                iter_v++;
                oldv2 = v;
                // Gradient step
                // TODO: special case when alpha_v = 0 -> S_v = I
                v = v + grad_v_step_size * (X.t()*u - S_v*v);

                // Proximal step
                v = (*prox_u)(v,prox_v_step_size);

                // Normalize with respect to S_v
                // Sometimes mat_norm(v,S_v) is so close to zero that v becomes NaN
                norm(v) > 0 ? v /= mat_norm(v, S_v) : v.zeros();

                in_v_tol = norm(v - oldv2) / norm(oldv2);
                // if(iter_v %100 == 0)
                MoMALogger::debug("---update v ") << iter_v << "---\n"
                                                  << "in_v_tol:" << in_v_tol
                                                  << "\t iter" << iter_v;
            }

            out_tol = norm(oldu1 - u) / norm(oldu1) + norm(oldv1 - v) / norm(oldv1);
            iter++;
            MoMALogger::debug("--Finish iter:") << iter << "---\n";
        }
    } else if (solver_type == Solver::FISTA){
        MoMALogger::error("FISTA is not provided yet!");
    } else {
        MoMALogger::error("Your choice of solver is not provided yet!");
    }

    MoMALogger::debug("==After the outer loop!==") << "out_tol:" << out_tol
                                                   << "\t iter" << iter;
}
