#include "moma_prox.h"

/*
* Prox base class
*/
NullProx::NullProx(){
    MoMALogger::debug("Initializing null proximal operator object");
}

arma::vec NullProx::operator()(const arma::vec &x, double l){
    return x;   // to be tested, return a reference might cause extra copying.
};

NullProx::~NullProx() {
    MoMALogger::debug("Releasing null proximal operator object");
};

/*
* Lasso
*/
Lasso::Lasso(){
    MoMALogger::debug("Initializing Lasso proximal operator object");
}

arma::vec Lasso::operator()(const arma::vec &x, double l){
    arma::vec absx = arma::abs(x);
    return arma::sign(x) % soft_thres_p(absx, l);
}

Lasso::~Lasso(){
    MoMALogger::debug("Releasing Lasso proximal operator object");
}

/*
* Non-negative Lasso
*/
NonNegativeLasso::NonNegativeLasso(){
    MoMALogger::debug("Initializing non-negative Lasso proximal operator object");
}

arma::vec NonNegativeLasso::operator()(const arma::vec &x, double l){
    return soft_thres_p(x,l);
}

NonNegativeLasso::~NonNegativeLasso(){
    MoMALogger::debug("Releasing non-negative Lasso proximal operator object");
}

/*
* SCAD
*/
SCAD::SCAD(double g){
    MoMALogger::debug("Initializing SCAD proximal operator object");
    if(g < 2){
        MoMALogger::error("Gamma for SCAD should be larger than 2!");
    };
    gamma = g;
}

SCAD::~SCAD(){
    MoMALogger::debug("Releasing SCAD proximal operator object");
}

arma::vec SCAD::operator()(const arma::vec &x, double l){
    int n = x.n_elem;
    double gl = gamma*l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    for (int i = 0; i < n; i++) // Probably need vectorization
    {
        // The implementation follows 
        // Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties
        // Jianqing Fan nd Runze Li
        // formula(2.8).
        z(i) = absx(i) > gl ? absx(i) : (absx(i) > 2 * l ? //(gamma-1)/(gamma-2) * THRES_P(absx(i),gamma*l/(gamma-1))
                                                ((gamma - 1) * absx(i) - gl)/ (gamma - 2)
                                                : THRES_P(absx(i),l)
                                                );
    }
    return z % sgnx;
}

arma::vec SCAD::vec_prox(const arma::vec &x, double l){
    int n = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    arma::umat D(x.n_elem,3,arma::fill::zeros);

    for(int i = 0; i < n; i++){
        arma::uword flag = absx(i) > gl ? 2 : (absx(i) > 2 * l ? 1: 0);
        D(i,flag) = 1;
    }

    // D.col(2) = absx > gl;
    // D.col(0) = absx <= 2 * l;
    // D.col(1) = arma::ones<arma::uvec>(n) - D.col(2) - D.col(0);
    z = D.col(0) % soft_thres_p(absx,l) + D.col(1) % ((gamma - 1) * absx - gl) / (gamma-2) + D.col(2) % absx;
    return sgnx % z;
}

/*
* Nonnegative SCAD
*/
NonNegativeSCAD::NonNegativeSCAD(double g):SCAD(g){
    MoMALogger::debug("Initializing non-negetive SCAD proximal operator object");
}

NonNegativeSCAD::~NonNegativeSCAD(){
    MoMALogger::debug("Releasing non-negetive SCAD proximal operator object");
}

arma::vec NonNegativeSCAD::operator()(const arma::vec &x, double l){
    int n = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    for (int i = 0; i < n; i++) // Probably need vectorization
    {
        // The implementation follows 
        // Variable Selection via Nonconcave Penalized Likelihood and its Oracle Properties
        // Jianqing Fan and Runze Li
        // formula(2.8).
        z(i) = x(i) > gl ? x(i): 
                                    (x(i) > 2 * l ? //(gamma-1)/(gamma-2) * THRES_P(absx(i),gamma*l/(gamma-1))
                                                ((gamma - 1) * x(i) - gl) / (gamma - 2)
                                                : THRES_P(x(i),l)
                                                );
    }
    return z;
}

/*
* MCP
*/
MCP::MCP(double g){
    MoMALogger::debug("Initializing MCP proximal operator object");
    if(g<1){
        MoMALogger::error("Gamma for MCP should be larger than 1!");
    }
    gamma=g;
}

MCP::~MCP(){
    MoMALogger::debug("Releasing MCP proximal operator object");
}

arma::vec MCP::operator()(const arma::vec &x, double l){
    int n = x.n_elem;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    for (int i = 0; i < n; i++)
    {
        // implementation follows lecture notes of Patrick Breheny
        // http://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-29.pdf
        // slide 19
        z(i) = absx(i) > gamma * l ? absx(i) : (gamma / (gamma - 1)) * THRES_P(absx(i),l);
    }
    return z % sgnx;    
}

arma::vec MCP::vec_prox(const arma::vec &x, double l){
    int n = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    arma::umat D(x.n_elem,2,arma::fill::zeros);
    for(int i = 0; i < n; i++){
        arma::uword flag = 0;
        if(absx(i) <= gl){
            flag = 1;
        }
        D(i,flag) = 1;
    }
    z = (gamma / (gamma - 1)) * (D.col(1) % soft_thres_p(absx,l)) + D.col(0) % absx;
    return sgnx % z;
}

/*
* Non-negative MCP
*/
NonNegativeMCP::NonNegativeMCP(double g):MCP(g){
    MoMALogger::debug("Initializing non-negative MCP proximal operator object");
}

NonNegativeMCP::~NonNegativeMCP(){
    MoMALogger::debug("Releasing non-negative MCP proximal operator object");
}

arma::vec NonNegativeMCP::operator()(const arma::vec &x, double l){
    int n = x.n_elem;
    arma::vec z(n);
    for (int i = 0; i < n; i++)
    {
        // implementation follows lecture notes of Patrick Breheny
        // http://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-29.pdf
        // slide 19
        z(i) = x(i) > gamma * l ? x(i) : (gamma / (gamma - 1)) * THRES_P(x(i),l);
    }
    return z;
}

/*
* Group lasso
*/
GrpLasso::GrpLasso(const arma::vec &grp):
        group(grp - arma::ones<arma::vec>(grp.n_elem)){   
        // takes in a factor `grp`, whose indices start with 1
    n_grp = grp.max();
    MoMALogger::debug("Initializing group lasso proximal operator object");
    D = arma::zeros<arma::umat>(n_grp,grp.n_elem);  // density will be 1/p = 1/x.n_elem
    for(int i = 0; i < grp.n_elem; i++){
        arma::uword g = group(i); // the i-th parameter is in g-th group. Note factor in R starts from 1
        D(g,i) = 1;
    }
}

GrpLasso::~GrpLasso(){
    MoMALogger::debug("Releasing non-negative group lasso proximal operator object");
}

arma::vec GrpLasso::operator()(const arma::vec &x, double l){
    // TODO: benchmark with simple looping!
    if(x.n_elem != group.n_elem){
        MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
    }
    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for(int i = 0; i < x.n_elem; i++){
        grp_norm(group(i)) += x(i)*x(i);
    }
    grp_norm = arma::sqrt(grp_norm);
    arma::vec grp_scale = soft_thres_p(grp_norm,l) / grp_norm;
    arma::vec scale(x.n_elem);
    for(int i = 0; i < x.n_elem; i++){
        scale(i) = grp_scale(group(i));
    }
    return x % scale;
}

arma::vec GrpLasso::vec_prox(const arma::vec &x, double l){
    arma::vec grp_norms = D.t() * arma::sqrt(D * arma::square(x)); // to_be_thres is of dimension p.
    return (x / grp_norms) % soft_thres_p(grp_norms,l);
};
/*
* Non-negative group lasso
*/
NonNegativeGrpLasso::NonNegativeGrpLasso(const arma::vec &grp):GrpLasso(grp){   // takes in a factor
    MoMALogger::debug("Initializing non-negative group lasso proximal operator object");
}

NonNegativeGrpLasso::~NonNegativeGrpLasso(){
    MoMALogger::debug("Releasing non-negative group lasso proximal operator object");
}

arma::vec NonNegativeGrpLasso::operator()(const arma::vec &x_, double l){
    // Reference: Proximal Methods for Hierarchical Sparse Coding, Lemma 11
    arma::vec x = soft_thres_p(x_,0);  // zero out negative entries
    if(x.n_elem != group.n_elem){
        MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
    }
    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for(int i = 0; i < x.n_elem; i++){
        grp_norm(group(i)) += x(i)*x(i);
    }
    grp_norm = arma::sqrt(grp_norm);
    arma::vec grp_scale = soft_thres_p(grp_norm,l) / grp_norm;
    arma::vec scale(x.n_elem);
    for(int i = 0; i < x.n_elem; i++){
        scale(i) = grp_scale(group(i));
    }
    return x % scale;
}

OrderedFusedLasso::OrderedFusedLasso(){
    MoMALogger::debug("Initializing a ordered fusion lasso proximal operator object");
}
OrderedFusedLasso::~OrderedFusedLasso(){
    MoMALogger::debug("Releasing a ordered fusion lasso proximal operator object");
}
arma::vec OrderedFusedLasso::operator()(const arma::vec &x, double l){
    FusedGroups fg(x);
    while(!fg.all_merged() && fg.next_lambda() < l){
        fg.merge();

    }
    return fg.find_beta_at(l);
}

/*
* Fusion lasso
*/
Fusion::Fusion(const arma::mat &input_w,bool input_ADMM,bool input_acc,double input_prox_eps){
    // w should be symmetric, and have zero diagonal elements.
    // We only store the upper half of the matrix w_ij, j >= i+1.
    // This wastes some space since the lower triangular part of weight is empty.
    prox_eps = input_prox_eps;
    ADMM = input_ADMM;
    acc = input_acc;
    int n_col = input_w.n_cols;
    int n_row = input_w.n_rows;
    if(n_col != n_row){
        MoMALogger::error("Weight matrix should have the same dimensions") << n_col << " and " << n_row;
    }
    weight.set_size(n_col,n_col);
    for(int i = 0; i < n_col; i++){
        for(int j = i + 1; j < n_col; j++){
            weight(i,j) = input_w(i,j);
        }
    }
    MoMALogger::debug("Initializing a fusion lasso proximal operator object");
};

Fusion::~Fusion(){
    MoMALogger::debug("Releasing a fusion lasso proximal operator object");
};

// Find the column sums and row sums of an upper triangular matrix,
// whose diagonal elements are all zero.
int tri_sums(const arma::mat &w, arma::vec &col_sums, arma::vec &row_sums, int n){
    // col_sums and row_sums should be initialzied by zeros.
    col_sums.zeros();
    row_sums.zeros();
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            row_sums(i) += w(i,j);
            col_sums(j) += w(i,j);
        }
    }
    return 0;
}

// This function sets lambda += (lambda - old_lambda) * step,
// and then set old_lambda = lambda.
int tri_momentum(arma::mat &lambda, arma::mat &old_lambda, double step, int n){
    arma::mat res(n,n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            double diff = lambda(i,j) - old_lambda(i,j);
            lambda(i,j) += step * diff;
            old_lambda(i,j) = lambda(i,j);
        }
    }
    return 0;
}

arma::vec Fusion::operator()(const arma::vec &x, double l){
    const int MAX_IT = 10000;
    arma::mat w = arma::trimatu(weight,1);
    int n = x.n_elem;
    if(n == 2){
        MoMALogger::error("Please use ordered fused lasso instead");
    }
    // AMA does not converge for current choice of nu.
    // Use ordianry ADMM instead.
    if(n < 20){
        ADMM = 1;
        acc = 0;    // TODO: when we provide acc ADMM
        MoMALogger::info("Use acc = 0, ADMM = 1 because AMA might not converge to optima.");
    }

    // beta subproblem: O(n)
    if(ADMM){
        // ADMM
        // Reference: Algorithm 5 in 
        // ADMM Algorithmic Regularization Paths for Sparse Statistical Machine Learning,
        // Yue Hu, Eric C. Chi and Genevera I. Allen
        // TODO: step size;i momentum

        // Using Genevera's paper notations
        const arma::vec &y = x;
        
        double y_bar = arma::mean(y);
        arma::mat z(n,n);
        arma::mat u(n,n);
        arma::vec b(n);
        arma::vec old_b;
        arma::mat old_z(n,n);
        arma::mat old_u(n,n);
        double old_alpha = 1;

        z.zeros();
        old_z.zeros();
        u.zeros();
        old_u.zeros();
        b.zeros();
        old_b.zeros();

        int cnt = 0;
        do{
            old_b = b;
            cnt++;
            arma::vec z_row_sums(n);
            arma::vec z_col_sums(n);
            arma::vec u_row_sums(n);
            arma::vec u_col_sums(n);
            tri_sums(u,u_col_sums,u_row_sums,n);
            tri_sums(z,z_col_sums,z_row_sums,n);
            // beta subproblem: O(n)
            for(int i = 0; i < n; i++){
                double part1 = z_row_sums(i) + u_row_sums(i);
                double part2 = z_col_sums(i) + u_col_sums(i);
                b(i) = ((y(i) + n * y_bar) + part1 - part2) / (n + 1);
            }
            // u and z subproblems: O(n(n-1)/2)
            for(int i = 0; i < n; i++){
                for(int j = i + 1; j < n; j++){
                    // z
                    double to_be_thres = b(i) - b(j) - u(i,j);
                    double scale = (1 - l * w(i,j) / std::abs(to_be_thres)); // TODO
                    if(scale < 0){
                        scale = 0;
                    };
                    z(i,j) = scale * to_be_thres;
                    // u
                    u(i,j) = u(i,j) + (z(i,j) - (b(i) - b(j)));
                }
            }
            if(acc){
                MoMALogger::error("Not provided yet");
                double alpha = (1 + std::sqrt(old_alpha)) / 2;
                z += (old_alpha / alpha) * (z - old_z);
                old_z = z;
                u += (old_alpha / alpha) * (u - old_u);
                old_u = u;
                // tri_momentum(z,old_z,old_alpha / alpha,n);
                // tri_momentum(u,old_u,old_alpha / alpha,n);
                old_alpha = alpha;
            }
        }while(arma::norm(old_b - b,2) / arma::norm(old_b,2) > prox_eps && cnt < MAX_IT);

        // Check convergence
        if(cnt == MAX_IT){
            MoMALogger::warning("No convergence in unordered fusion lasso prox (ADMM).");
        }else{
            MoMALogger::debug("ADMM converges.");
        }
        // TODO: shrink stepsize, as is done in the paper
        return b;
    }
    else{
        // AMA
        // Reference: Algorithm 3 Fast AMA in
        // Splitting Methods for Convex Clustering
        // Eric C. Chi and Kenneth Lange†
        int n = x.n_elem;
        // Choosing nu is not trivial. See Proposition 4.1 in 
        // Splitting Methods for Convex Clustering, Chi and Range

        // Find the degrees of nodes
        arma::vec deg(n);
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                if(w(i,j) > 0){
                    deg(i)++;
                    deg(j)++;
                }
            }
        }
        // Find the degree of edges, which are the sums of their nodes.
        int max_edge_deg = -1;
        for(int i = 0; i < n; i++){
            for(int j = i + 1; j < n; j++){
                if(w(i,j) > 0){
                    max_edge_deg = std::max(deg(i)+deg(j),(double)max_edge_deg);
                }
            }
        }
        double nu = 1.0 / (std::min(n,max_edge_deg));
        arma::vec u(n);
        arma::mat lambda(n,n);

        // Initialze
        u.zeros();
        lambda.zeros();
        double old_alpha = 1;
        arma::mat old_lambda(n,n);
        old_lambda.zeros();
        arma::vec old_u;
        int cnt = 0;

        // Start iterating
        do{
            cnt++;
            // lambda subproblem
            for(int i = 0; i < n; i++){
                for(int j = i + 1; j < n; j++){
                    lambda(i,j) = lambda(i,j) - nu * (u(i) - u(j));
                    // project onto the interval [ -w_ij*lambda_ij, w_ij*lambda_ij ]
                    if(std::abs(lambda(i,j)) > l * w(i,j)){
                        lambda(i,j) = l * w(i,j) * lambda(i,j) / std::abs(lambda(i,j));
                    }
                }
            }

            // u subproblem
            old_u = u;
            arma::vec lambda_row_sums(n);
            arma::vec lambda_col_sums(n);
            tri_sums(lambda,lambda_col_sums,lambda_row_sums,n);
            for(int i = 0; i < n; i++){
                double part1 = lambda_row_sums(i);
                double part2 = lambda_col_sums(i);
                u(i) = x(i) + part1 - part2;
            }
            if(acc){// Momemtum step
                double alpha = (1 + std::sqrt(old_alpha)) / 2;
                tri_momentum(lambda,old_lambda,old_alpha / alpha,n);
                old_alpha = alpha;
            }
        }while(arma::norm(u-old_u,2) / arma::norm(old_u,2) > prox_eps && cnt < MAX_IT);

        // Check convergence
        if(cnt == MAX_IT){
            MoMALogger::warning("No convergence in unordered fusion lasso prox (AMA).");
        }else{
            MoMALogger::debug("AMA converges.");
        }
        return u;
    }
}
