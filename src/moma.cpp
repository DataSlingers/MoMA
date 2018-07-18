// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "moma.h"

enum class Solver{
    ISTA,
    FISTA,
    APP_ISTA
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
    double alpha_u;
    double alpha_v;
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
        double i_alpha_u,    // Smoothing levels
        double i_alpha_v,
        /*
         * Algorithm parameters:
         */
        double i_EPS,
        arma::uword i_MAX_ITER,
        std::string i_solver):alpha_u(i_alpha_u),alpha_v(i_alpha_v),X(X_) // const reference must be passed to initializer list
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
        S_u += alpha_u * Omega_u;
        S_v += alpha_v * Omega_v;
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
        if(mat_norm(u,S_u) != 0){  // Normalize one more time just in case
            u = u / mat_norm(u,S_u);
        }
        if(mat_norm(v,S_v) != 0){
            v = v / mat_norm(v,S_v);
        }
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
    Prox* res = new NullProx();
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
    else if(s.compare("ORDEREDFUSED") == 0){
        if(nonneg)
            MoMALogger::error("Non-negative ordered fused lasso is not implemented!");
        else
            res = new OrderedFusion();
    }
    else
        MoMALogger::warning("Your sparse penalty is not provided by us/specified by you! Use `NullProx` by default");
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
        MoMALogger::info("Model info=========")
            <<"n:" <<n
            <<"p:" << p;
        MoMALogger::info("Start fitting.");

        // store the value of u and v at the start of outer loop, hence call it oldu1
        arma::vec oldu1 = arma::zeros<arma::vec>(n);
        arma::vec oldv1 = arma::zeros<arma::vec>(p);
        // store the value of u and v at the start of inner loop
        arma::vec oldu2 = arma::zeros<arma::vec>(n);
        arma::vec oldv2 = arma::zeros<arma::vec>(p);

        // number of iteration
        int iter = 0;
        int iter_u = 0;
        int iter_v = 0;

        // stopping tolerance
        double in_u_tol = 1;   // tolerance for inner loop of u updates
        double in_v_tol = 1;   // tolerance for inner loop of v updates
        double out_tol = 1;    // that of outer loop

        if (solver_type == Solver::FISTA){
            MoMALogger::info("Running FISTA!");
            MoMALogger::debug("==Before the loop: training setup==")
                    << "\titer" << iter
                    << "\tEPS:" << EPS
                    << "\tMAX_ITER:" << MAX_ITER;
            while (out_tol > EPS && iter < MAX_ITER)
            {
                oldu1 = u;
                oldv1 = v;
                in_u_tol = 1;
                in_v_tol = 1;
                iter_u = 0;
                iter_v = 0;

                /***********
                * Update of u
                ************/
                double t = 1;   // momemtum stepsize
                double mn = 0;  // matrix norm
                while (in_u_tol > EPS && iter_u < MAX_ITER)
                {
                    iter_u++;
                    oldu2 = u;
                    double oldt = t;
                    t = 0.5 * (1 + sqrt(1 + 4 * oldt*oldt));
                    // gradient step
                    if(alpha_u == 0.0){
                        u = u + grad_u_step_size * (X*v - u);
                    }else{
                        u = u + grad_u_step_size * (X*v - S_u*u);
                    }
                    // proxiaml step
                    u = (*prox_u)(u,prox_u_step_size);
                    // momemtum step
                    u = u + (oldt - 1) / t * (u - oldu2);
                    // find torlerance
                    in_u_tol = norm(u - oldu2) / norm(oldu2);
                    MoMALogger::debug("u ") << iter_v << "---" << "% of change " << in_u_tol;
                }
                // nomalize w.r.t S_u
                mn = mat_norm(u, S_u);
                mn > 0 ? u /= mn : u.zeros();
                MoMALogger::debug("mat_norm is ")  << mn;


                /***********
                * Update of v
                ************/
                // restore
                t = 1;
                while (in_v_tol > EPS && iter_v < MAX_ITER)
                {
                    iter_v++;
                    oldv2 = v;
                    double oldt = t;
                    t = 0.5 * (1 + sqrt(1 + 4 * oldt*oldt));
                    // gradient step
                    if(alpha_v == 0.0){
                        v = v + grad_v_step_size * (X.t()*u - v);
                    }else{
                       v = v + grad_v_step_size * (X.t()*u - S_v*v);
                    }
                    // proximal step
                    v = (*prox_v)(v,prox_v_step_size);
                    // momemtum step
                    v = v + (oldt - 1) / t * (v - oldv2);
                    // find tolerance
                    in_v_tol = norm(v - oldv2) / norm(oldv2);
                    MoMALogger::debug("v ") << iter_v << "---" << "% of change " << in_v_tol;
                }
                // normalize w.r.t. S_v
                mn = mat_norm(v, S_v);
                mn > 0 ? v /= mn : v.zeros();
                MoMALogger::debug("mat_norm is ") << mn;

                // Output info
                out_tol = norm(oldu1 - u) / norm(oldu1) + norm(oldv1 - v) / norm(oldv1);
                iter++;
                MoMALogger::info("--Finish iter:") << iter << "---";
            }
        }
        else if (solver_type == Solver::ISTA) {
            MoMALogger::info("Running ISTA!");
            MoMALogger::debug("==Before the loop: training setup==")
                    << "\titer" << iter
                    << "\tEPS:" << EPS
                    << "\tMAX_ITER:" << MAX_ITER;
            while (out_tol > EPS && iter < MAX_ITER)
            {
                oldu1 = u;
                oldv1 = v;
                in_u_tol = 1;
                in_v_tol = 1;
                iter_u = 0;
                iter_v = 0;

                /***********
                * Update of u
                ************/
                double mn = 0;
                while (in_u_tol > EPS && iter_u < MAX_ITER)
                {
                    iter_u++;
                    oldu2 = u;
                    // gradient step
                    if(alpha_u == 0.0){
                        u = u + grad_u_step_size * (X*v - u);
                    }else{
                        u = u + grad_u_step_size * (X*v - S_u*u);
                    }
                    // proxiaml step
                    u = (*prox_u)(u,prox_u_step_size);
                    // find tolerance
                    in_u_tol = norm(u - oldu2) / norm(oldu2);
                    MoMALogger::debug("u ") << iter_u << "--"<< "% of change " << in_u_tol;
                }
                // nomalize w.r.t S_u
                mn = mat_norm(u, S_u);
                mn > 0 ? u /= mn : u.zeros();
                MoMALogger::debug("mat_norm is ") << mn;


                /***********
                * Update of v
                ************/
                while (in_v_tol > EPS && iter_v < MAX_ITER)
                {
                    iter_v++;
                    oldv2 = v;
                    // gradient step
                    if(alpha_v == 0.0){
                        v = v + grad_v_step_size * (X.t()*u - v);
                    }else{
                       v = v + grad_v_step_size * (X.t()*u - S_v*v);
                    }
                    // proximal step
                    v = (*prox_v)(v,prox_v_step_size);
                    // find tolerance
                    in_v_tol = norm(v - oldv2) / norm(oldv2);
                    MoMALogger::debug("v ") << iter_v << "---" << "% of change " << in_v_tol;
                }
                // nomalize w.r.t S_v
                mn = mat_norm(v, S_v);
                mn > 0 ? v /= mn : v.zeros();
                MoMALogger::debug("mat_norm is ") << mn;

                out_tol = norm(oldu1 - u) / norm(oldu1) + norm(oldv1 - v) / norm(oldv1);
                iter++;
                MoMALogger::info("--Finish iter:") << iter << "---";
            }
        }
        else{
            MoMALogger::error("Your choice of solver is not provided yet!");
        }
        MoMALogger::debug("==After the outer loop!==") 
                   << "out_tol:" << out_tol << "\t iter" << iter;
        if(iter == MAX_ITER)
            MoMALogger::warning("No convergence!");
}
