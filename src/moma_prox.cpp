#include "moma_prox.h"

/*
* Prox base class
*/
Prox::Prox(){
    MoMALogger::debug("Initializing proximal operator object");
}

arma::vec Prox::operator()(const arma::vec &x, double l){
    return x;   // to be tested, return a reference might cause extra copying.
};

Prox::~Prox() {
    MoMALogger::debug("Releasing proximal operator object");
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
    if(g < 2)
        MoMALogger::error("Gamma for SCAD should be larger than 2!");
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
    if(g<1) 
        MoMALogger::error("Gamma for MCP should be larger than 1!");
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
        if(absx(i) <= gl)
            flag = 1;   
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
GrpLasso::GrpLasso(const arma::vec &grp):group(grp - arma::ones<arma::vec>(grp.n_elem)){   // takes in a factor grp, whose indeces start with 1
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
    if(x.n_elem != group.n_elem)
        MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
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
    if(x.n_elem != group.n_elem)
    MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
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
