#include "moma_prox.h"

/*
* L1 Trend Filtering
*/

// A second difference matrix
// [0 0 ... 1 -2 1 ... 0 0]
arma::mat sec_diff_mat(int m){
    if(m < 3){
        MoMALogger::error("A second difference matrix should have more that 3 columns.");
    }
    arma::mat D = arma::zeros<arma::mat>(m-2,m);
    for(int i = 0; i < m-2; i++){
        D(i,i) = 1;
        D(i,i+1) = -2;
        D(i,i+2) = 1;
    }
    return D;
}

L1TrendFiltering::L1TrendFiltering(){
    MoMALogger::debug("Initializing a L1 trend filtering proximal operator object");
}

L1TrendFiltering::~L1TrendFiltering(){
    MoMALogger::debug("Releasing a L1 trend filtering proximal operator object");
}

// Find the sum of dual residual and centering residual
double res(
    const arma::mat &DDt, const arma::vec &Dy,
    const arma::vec &mu1,const arma::vec &mu2,
    const arma::vec &f1,const arma::vec &f2,double t,
    const arma::vec &nu){
    arma::vec res_cent = - (mu1 % f1) - (mu2 % f2) - 1/t;
    // NOTE: this approximates formulea 16 in the paper
    arma::vec res_dual = DDt * nu - Dy + mu1 - mu2;
    return arma::as_scalar(arma::norm(res_cent) + arma::norm(res_dual));
}

double init_stepsize(const arma::vec &mu,const arma::vec &dmu,double step){
    if(step <= 0){
        MoMALogger::error("step should > 0 in init_stepsize.");
    }

    arma::uvec idx = (dmu < 0);

    if(arma::sum(idx) > 0){
        // the largest step size that does not negate any elements of mu
        // i.e., max(step) s.t. mu + step * dmu > 0 for all elements
        arma::vec prop = -(idx % (mu / dmu));

        double max_legit_ss = step;
        for(int i = 0; i < mu.n_elem; i++){
            if(prop(i) > 0 && prop(i) < max_legit_ss){
                max_legit_ss = 0.99 * prop(i);
            }
        }
        return max_legit_ss;
    }
    else{
        return step;
    }
}

arma::vec L1TrendFiltering::operator()(const arma::vec &x, double l){
    int MAX_ITER = 20;
    int MAX_BT_ITER = 5;
    double prox_eps = 1e-10;

    // The backtracking parameters
    // shrink stepsize by `bata`
    // if f(x + stepsize * dx) >= (1 - alpha * step) * f(x)
    // Ref: http://www.stat.cmu.edu/~ryantibs/convexopt-F15/lectures/16-primal-dual.pdf page 12
    double alpha = 0.01;
    double beta = 0.5;

    int n = x.n_elem;
    const arma::vec &y = x; 

    // Commonly used mat and vec
    arma::mat D = sec_diff_mat(n);
    arma::mat DDt = D * D.t();
    arma::vec Dy = D * y;

    // Step size
    double step;

    // Initialzation with a strictly feasible point
    // Ref: l1 Trend Filtering by Seung-Jean Kim, Kwangmoo Koh
    // Stephen Boyd and Dimitry Gorinevsky
    arma::vec nu = arma::zeros<arma::vec> (n-2);        // The dual variable, y - D^t*nu is what we need
    arma::vec mu1 = arma::ones<arma::vec> (n-2);        // Multipliers to solve the dual problem
    arma::vec mu2 = arma::ones<arma::vec> (n-2);

    arma::vec new_nu;
    arma::vec new_mu1;
    arma::vec new_mu2;

    // The inequality constraints, all elements should be non-positive
    arma::vec f1 = nu - l;
    arma::vec f2 = -nu - l;

    arma::vec new_f1;
    arma::vec new_f2;

    // We solve the dual problems with increasing t
    double t = 1e-10;

    // Surrogate duality gap = - sum_i {u_i * h_i(x)}
    double gap;

    int iter = 0;
    for(; iter < MAX_ITER; iter++){

        // Surrogate duality gap
        // Ref: PPT page 10
        gap = - (sum(mu1 % f1 + mu2 % f2));
        if(gap < prox_eps){
            break;
        }

        // Ref: PPT page 11
        t = 4 * (n-2) / gap;

        arma::mat part1 = DDt - arma::diagmat(mu1/f1+mu2/f2);       // A band matrix
        arma::vec part2 = -DDt*nu + Dy + (1/t) / f1 - (1/t) / f2;

        // Update directions
        arma::vec dnu = arma::solve(part1,part2);                   // RcppArmadillo optimizes with band matrix automatically
        arma::vec dmu1 = -(mu1 + ((1/t) + dnu % mu1) / f1);
        arma::vec dmu2 = -(mu2 + ((1/t) + dnu % mu2) / f2);

        // Choose step size by back tracking
        // Initialize with the largest stepsize not making mu's negative
        step = init_stepsize(mu1,dmu1,1);
        step = init_stepsize(mu2,dmu2,step);

        for(int iter_bt = 0; iter_bt < MAX_BT_ITER; iter_bt++){

            new_nu = nu + step * dnu;
            new_mu1 = mu1 + step * dmu1;
            new_mu2 = mu2 + step * dmu2;

            // new inequality constraints
            new_f1 = new_nu - l;
            new_f2 = -new_nu - l;

            // check step size, exit if all of the conditions are met
            // 1. f1<0 and f2<0
            // 2. loss function should be less than a reference linear function
            if(arma::sum(new_f1 > 0) > 0 || arma::sum(new_f2 > 0) > 0 ||
                res(DDt,Dy,new_mu1,new_mu2,new_f1,new_f2,t,new_nu) >(1 - alpha * step) * res(DDt,Dy,mu1,mu2,f1,f2,t,nu)){
                        step *= beta;
            }
            else{
                break;
            }
        }

        // Update variable with selected stepsize
        nu = new_nu;
        mu1 = new_mu1;
        mu2 = new_mu2;
        f1 = new_f1;
        f2 = new_f2;
    }

    if(iter == MAX_ITER){
        MoMALogger::warning("No convergence in L1 trend filtering solver. Surrogate duality gap = ") 
                    << gap 
                    << ".";
    }

    return y - D.t() * nu;
}
