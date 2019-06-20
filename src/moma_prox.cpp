#include "moma_prox.h"

/*
 * Prox base class
 */
NullProx::NullProx()
{
    MoMALogger::debug("Initializing null proximal operator object");
}

arma::vec NullProx::operator()(const arma::vec &x, double l)
{
    return x;  // TODO: to be tested, return a reference might cause extra
               // copying.
};

NullProx::~NullProx()
{
    MoMALogger::debug("Releasing null proximal operator object");
};

int NullProx::df(const arma::vec &x)
{
    return x.n_elem;
}

/*
 * Lasso
 */
Lasso::Lasso()
{
    MoMALogger::debug("Initializing Lasso proximal operator object");
}

arma::vec Lasso::operator()(const arma::vec &x, double l)
{
    arma::vec absx = arma::abs(x);
    return arma::sign(x) % soft_thres_p(absx, l);
}

Lasso::~Lasso()
{
    MoMALogger::debug("Releasing Lasso proximal operator object");
}

int Lasso::df(const arma::vec &x)
{
    return arma::sum(x != 0.0);
}

/*
 * SLOPE - Sorted L-One Penalized Estimation
 */
SLOPE::SLOPE(int dim)
{
    // BH-type rule
    // lambda_BH(i) = Phi_inv ()
    lambda.resize(dim);
    double q = 0.05;
    for (int i = 0; i < dim; i++)
    {
        lambda(i) = R::qnorm(1 - (i + 1) * q / (2 * dim), 0.0, 1.0, 1, 0);
    }

    MoMALogger::debug("Initializing SLOPE proximal operator object");
}

arma::vec SLOPE::operator()(const arma::vec &x, double l)
{
    int n           = x.n_elem;
    arma::vec x_sgn = arma::sign(x);
    arma::vec x_abs = arma::abs(x);

    arma::uvec order = arma::sort_index(x_abs, "descend");

    arma::vec ordered_absx(n);
    //  arma::vec ordered_abs(n);
    for (int i = 0; i < n; i++)
    {
        ordered_absx(i) = x_abs(order(i));
    }

    arma::vec scratch(x.n_elem);
    evaluateProx(ordered_absx, lambda * l, scratch, x.n_elem, order);
    return scratch % x_sgn;
}

SLOPE::~SLOPE()
{
    MoMALogger::debug("Releasing SLOPE proximal operator object");
}

int SLOPE::df(const arma::vec &x)
{
    return arma::sum(x != 0.0);
}

/*
 * Non-negative Lasso
 */
NonNegativeLasso::NonNegativeLasso()
{
    MoMALogger::debug("Initializing non-negative Lasso proximal operator object");
}

arma::vec NonNegativeLasso::operator()(const arma::vec &x, double l)
{
    return soft_thres_p(x, l);
}

NonNegativeLasso::~NonNegativeLasso()
{
    MoMALogger::debug("Releasing non-negative Lasso proximal operator object");
}

int NonNegativeLasso::df(const arma::vec &x)
{
    return arma::sum(x != 0.0);
}

/*
 * SCAD
 */
SCAD::SCAD(double g)
{
    MoMALogger::debug("Initializing SCAD proximal operator object");
    if (g <= 2)
    {
        MoMALogger::error("Gamma for SCAD should be larger than 2!");
    };
    gamma = g;
}

SCAD::~SCAD()
{
    MoMALogger::debug("Releasing SCAD proximal operator object");
}

arma::vec SCAD::operator()(const arma::vec &x, double l)
{
    int n     = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    for (int i = 0; i < n; i++)  // Probably need vectorization
    {
        // The implementation follows
        // Variable Selection via Nonconcave Penalized Likelihood and its Oracle
        // Properties Jianqing Fan nd Runze Li formula(2.8).
        z(i) = absx(i) > gl
                   ? absx(i)
                   : (absx(i) > 2 * l ?  //(gamma-1)/(gamma-2) * THRES_P(absx(i),gamma*l/(gamma-1))
                          ((gamma - 1) * absx(i) - gl) / (gamma - 2)
                                      : THRES_P(absx(i), l));
    }
    return z % sgnx;
}

arma::vec SCAD::vec_prox(const arma::vec &x, double l)
{
    int n     = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    arma::umat D(x.n_elem, 3, arma::fill::zeros);

    for (int i = 0; i < n; i++)
    {
        arma::uword flag = absx(i) > gl ? 2 : (absx(i) > 2 * l ? 1 : 0);
        D(i, flag)       = 1;
    }

    // D.col(2) = absx > gl;
    // D.col(0) = absx <= 2 * l;
    // D.col(1) = arma::ones<arma::uvec>(n) - D.col(2) - D.col(0);
    z = D.col(0) % soft_thres_p(absx, l) + D.col(1) % ((gamma - 1) * absx - gl) / (gamma - 2) +
        D.col(2) % absx;
    return sgnx % z;
}

int SCAD::df(const arma::vec &x)
{
    // An approximation
    return arma::sum(x != 0.0);
}

/*
 * Nonnegative SCAD
 */
NonNegativeSCAD::NonNegativeSCAD(double g) : SCAD(g)
{
    MoMALogger::debug("Initializing non-negative SCAD proximal operator object");
}

NonNegativeSCAD::~NonNegativeSCAD()
{
    MoMALogger::debug("Releasing non-negative SCAD proximal operator object");
}

arma::vec NonNegativeSCAD::operator()(const arma::vec &x, double l)
{
    int n     = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    for (int i = 0; i < n; i++)  // Probably need vectorization
    {
        // The implementation follows
        // Variable Selection via Nonconcave Penalized Likelihood and its Oracle
        // Properties Jianqing Fan and Runze Li formula(2.8).
        z(i) = x(i) > gl
                   ? x(i)
                   : (x(i) > 2 * l ?  //(gamma-1)/(gamma-2) * THRES_P(absx(i),gamma*l/(gamma-1))
                          ((gamma - 1) * x(i) - gl) / (gamma - 2)
                                   : THRES_P(x(i), l));
    }
    return z;
}

int NonNegativeSCAD::df(const arma::vec &x)
{
    // An approximation
    return arma::sum(x != 0.0);
}

/*
 * MCP
 */
MCP::MCP(double g)
{
    MoMALogger::debug("Initializing MCP proximal operator object");
    if (g <= 1)
    {
        MoMALogger::error("Gamma for MCP should be larger than 1!");
    }
    gamma = g;
}

MCP::~MCP()
{
    MoMALogger::debug("Releasing MCP proximal operator object");
}

arma::vec MCP::operator()(const arma::vec &x, double l)
{
    int n = x.n_elem;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    for (int i = 0; i < n; i++)
    {
        // implementation follows lecture notes of Patrick Breheny
        // http://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-29.pdf
        // slide 19
        z(i) = absx(i) > gamma * l ? absx(i) : (gamma / (gamma - 1)) * THRES_P(absx(i), l);
    }
    return z % sgnx;
}

arma::vec MCP::vec_prox(const arma::vec &x, double l)
{
    int n     = x.n_elem;
    double gl = gamma * l;
    arma::vec z(n);
    arma::vec absx = arma::abs(x);
    arma::vec sgnx = arma::sign(x);
    arma::umat D(x.n_elem, 2, arma::fill::zeros);
    for (int i = 0; i < n; i++)
    {
        arma::uword flag = 0;
        if (absx(i) <= gl)
        {
            flag = 1;
        }
        D(i, flag) = 1;
    }
    z = (gamma / (gamma - 1)) * (D.col(1) % soft_thres_p(absx, l)) + D.col(0) % absx;
    return sgnx % z;
}

int MCP::df(const arma::vec &x)
{
    // An approximation
    return arma::sum(x != 0.0);
}

/*
 * Non-negative MCP
 */
NonNegativeMCP::NonNegativeMCP(double g) : MCP(g)
{
    MoMALogger::debug("Initializing non-negative MCP proximal operator object");
}

NonNegativeMCP::~NonNegativeMCP()
{
    MoMALogger::debug("Releasing non-negative MCP proximal operator object");
}

arma::vec NonNegativeMCP::operator()(const arma::vec &x, double l)
{
    int n = x.n_elem;
    arma::vec z(n);
    for (int i = 0; i < n; i++)
    {
        // implementation follows lecture notes of Patrick Breheny
        // http://myweb.uiowa.edu/pbreheny/7600/s16/notes/2-29.pdf
        // slide 19
        z(i) = x(i) > gamma * l ? x(i) : (gamma / (gamma - 1)) * THRES_P(x(i), l);
    }
    return z;
}

int NonNegativeMCP::df(const arma::vec &x)
{
    // An approximation
    return arma::sum(x != 0.0);
}

/*
 * Group lasso
 */
GrpLasso::GrpLasso(const arma::vec &grp) : group(grp - arma::ones<arma::vec>(grp.n_elem))
{
    // takes in a factor `grp`, whose indices start with 1
    n_grp = grp.max();
    MoMALogger::debug("Initializing group lasso proximal operator object");
}

GrpLasso::~GrpLasso()
{
    MoMALogger::debug("Releasing non-negative group lasso proximal operator object");
}

arma::vec GrpLasso::operator()(const arma::vec &x, double l)
{
    // TODO: benchmark with simple looping!
    if (x.n_elem != group.n_elem)
    {
        MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
    }
    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for (int i = 0; i < x.n_elem; i++)
    {
        grp_norm(group(i)) += x(i) * x(i);
    }
    grp_norm            = arma::sqrt(grp_norm);
    arma::vec grp_scale = soft_thres_p(grp_norm, l) / grp_norm;
    arma::vec scale(x.n_elem);
    for (int i = 0; i < x.n_elem; i++)
    {
        scale(i) = grp_scale(group(i));
    }
    return x % scale;
}

int GrpLasso::df(const arma::vec &x)
{
    // Ref: Equation (6.3) of
    // Yuan, Ming, and Yi Lin.
    // "Model selection and estimation in regression with grouped variables."
    // Journal of the Royal Statistical Society: Series B (Statistical
    // Methodology) 68.1 (2006): 49-67.

    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for (int i = 0; i < x.n_elem; i++)
    {
        grp_norm(group(i)) += x(i) * x(i);
    }
    // Approxiamte df = number of groups that are not zero + num_para - num_group
    // The unbiased estimate is a bit complicated, involving finding OLS
    // estimates.
    return arma::sum(grp_norm != 0.0) + x.n_elem - n_grp;
}

/*
 * Non-negative group lasso
 */
NonNegativeGrpLasso::NonNegativeGrpLasso(const arma::vec &grp) : GrpLasso(grp)
{  // takes in a factor
    MoMALogger::debug("Initializing non-negative group lasso proximal operator object");
}

NonNegativeGrpLasso::~NonNegativeGrpLasso()
{
    MoMALogger::debug("Releasing non-negative group lasso proximal operator object");
}

arma::vec NonNegativeGrpLasso::operator()(const arma::vec &x_, double l)
{
    // Reference: Proximal Methods for Hierarchical Sparse Coding, Lemma 11
    arma::vec x = soft_thres_p(x_, 0);  // zero out negative entries
    if (x.n_elem != group.n_elem)
    {
        MoMALogger::debug("Wrong dimension: x dim is") << x.n_elem << "but we take" << group.n_elem;
    }
    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for (int i = 0; i < x.n_elem; i++)
    {
        grp_norm(group(i)) += x(i) * x(i);
    }
    grp_norm            = arma::sqrt(grp_norm);
    arma::vec grp_scale = soft_thres_p(grp_norm, l) / grp_norm;
    arma::vec scale(x.n_elem);
    for (int i = 0; i < x.n_elem; i++)
    {
        scale(i) = grp_scale(group(i));
    }
    return x % scale;
}

int NonNegativeGrpLasso::df(const arma::vec &x)
{
    arma::vec grp_norm = arma::zeros<arma::vec>(n_grp);
    for (int i = 0; i < x.n_elem; i++)
    {
        grp_norm(group(i)) += x(i) * x(i);
    }
    // Approxiamte df = number of groups that are not zero + num_para - num_group
    // The unbiased estimate is a bit complicated, involving find OLS estimates.
    return arma::sum(grp_norm != 0.0) + x.n_elem - n_grp;
}

/*
 * Ordered fused lasso
 */
OrderedFusedLasso::OrderedFusedLasso()
{
    MoMALogger::debug("Initializing a ordered fusion lasso proximal operator object");
}

OrderedFusedLasso::~OrderedFusedLasso()
{
    MoMALogger::debug("Releasing a ordered fusion lasso proximal operator object");
}

arma::vec OrderedFusedLasso::operator()(const arma::vec &x, double l)
{
    FusedGroups fg(x);
    while (!fg.all_merged() && fg.next_lambda() < l)
    {
        fg.merge();
    }
    return fg.find_beta_at(l);
}

int OrderedFusedLasso::df(const arma::vec &x)
{
    // Ref:
    // Table 2 of Tibshirani, Robert.
    // "Regression shrinkage and selection via the lasso: a retrospective."
    // Journal of the Royal Statistical Society: Series B (Statistical
    // Methodology) 73.3 (2011): 273-282.
    if (x.n_elem < 2)
    {
        MoMALogger::error("x must not be a scalar");
    }
    int df = 1;
    // Cound the number of transitions, which is the number of fused groups
    for (int i = 1; i < x.n_elem; i++)
    {
        if (abs(x(i) - x(i - 1)) > 1e-10)
        {
            df++;
        }
    }
    return df;
}

/*
 * Ordered fused lasso-dynamic programming approach
 */
OrderedFusedLassoDP::OrderedFusedLassoDP()
{
    MoMALogger::debug("Initializing a ordered fusion lasso proximal operator object (DP)");
}

OrderedFusedLassoDP::~OrderedFusedLassoDP()
{
    MoMALogger::debug("Releasing a ordered fusion lasso proximal operator object (DP)");
}

arma::vec OrderedFusedLassoDP::operator()(const arma::vec &x, double l)
{
    return myflsadp(x, l, MOMA_FUSEDLASSODP_BUFFERSIZE);
}

/*
 * Sparse fused lasso
 */
SparseFusedLasso::SparseFusedLasso(double i_lambda2) : lambda2(i_lambda2)
{
    MoMALogger::debug("Initializing a sparse fused lasso proximal operator object");
}

SparseFusedLasso::~SparseFusedLasso()
{
    MoMALogger::debug("Releasing a sparse fused lasso proximal operator object");
}

arma::vec SparseFusedLasso::operator()(const arma::vec &x, double l)
{
    arma::vec tmp = fg(x, l);
    return soft_thres(tmp, lambda2);
}

int SparseFusedLasso::df(const arma::vec &x)
{
    // Ref:
    // Table 2 of Tibshirani, Robert.
    // "Regression shrinkage and selection via the lasso: a retrospective."
    // Journal of the Royal Statistical Society: Series B (Statistical
    // Methodology) 73.3 (2011): 273-282.

    if (x.n_elem < 2)
    {
        MoMALogger::error("x must not be a scalar");
    }
    int df = x(0) != 0;
    // the number of non-zero fused groups
    for (int i = 1; i < x.n_elem; i++)
    {
        if (x(i) != 0 && abs(x(i) - x(i - 1)) > 1e-10)
        {
            df++;
        }
    }
    return df;
}

/*
 * Fusion lasso
 */

// From matrix index to vector index of an upper triangular matrix
int tri_idx(int i, int j, int n)
{
    return (2 * n - i - 1) * i / 2 + j - i - 1;
}

Fusion::Fusion(const arma::mat &input_w, bool input_ADMM, bool input_acc, double input_prox_eps)
{
    // w should be symmetric, and have zero diagonal elements.
    // We only store the upper half of the matrix w_ij, j >= i+1.
    // This wastes some space since the lower triangular part of weight is empty.
    prox_eps  = input_prox_eps;
    ADMM      = input_ADMM;
    acc       = input_acc;
    int n_col = input_w.n_cols;
    int n_row = input_w.n_rows;
    if (n_col != n_row)
    {
        MoMALogger::error("Weight matrix should be square: ") << n_col << " and " << n_row;
    }
    start_point.set_size(n_col);
    start_point.zeros();
    weight.set_size(n_col * (n_col - 1) / 2);

    for (int i = 0; i < n_col; i++)
    {
        for (int j = i + 1; j < n_col; j++)
        {
            int k     = tri_idx(i, j, n_col);
            weight(k) = input_w(i, j);
        }
    }
    MoMALogger::debug("Initializing a fusion lasso proximal operator object");
};

Fusion::~Fusion()
{
    MoMALogger::debug("Releasing a fusion lasso proximal operator object");
};

// Find the column sums and row sums of an upper triangular matrix,
// whose diagonal elements are all zero.
int tri_sums(const arma::vec &w, arma::vec &col_sums, arma::vec &row_sums, int n)
{
    // col_sums and row_sums should be initialzied by zeros.
    col_sums.zeros();
    row_sums.zeros();
    int cnt = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            row_sums(i) += w(cnt);
            col_sums(j) += w(cnt);
            cnt++;
        }
    }
    return 0;
}

int Fusion::df(const arma::vec &x)
{
    // First construct a graph, where any nodes are connected if
    // they have the same value. Then a good intuition is to find
    // the number of connected components, which is the multiplicity
    // of 0 as an eigenvalue of the laplacian matrix of the graph. Way too
    // complicated.

    // There could be cases when two groups happen to have the
    // same value, but are not actually "fused".

    // Let's just count the number of different values.
    arma::vec uq = arma::unique(x);
    return arma::sum(uq > 0);
}

// This function sets lambda += (lambda - old_lambda) * step,
// and then set old_lambda = lambda.
int tri_momentum(arma::mat &lambda, arma::mat &old_lambda, double step, int n)
{
    lambda += step * (lambda - old_lambda);
    old_lambda = lambda;
    return 0;
}

arma::vec Fusion::operator()(const arma::vec &x, double l)
{
    const int MAX_IT = 10000;
    arma::vec &w     = weight;
    int n            = x.n_elem;
    if (n == 2)
    {
        MoMALogger::error("Please use ordered fused lasso instead");
    }

    // beta subproblem: O(n)
    if (ADMM)
    {
        // ADMM
        // Reference: Algorithm 5 in
        // ADMM Algorithmic Regularization Paths for Sparse Statistical Machine
        // Learning, Yue Hu, Eric C. Chi and Genevera I. Allen

        // Using Genevera's paper notations
        MoMALogger::debug("Running ADMM.");
        const arma::vec &y = x;

        double y_bar = arma::mean(y);
        arma::vec z(n * (n - 1) / 2);
        arma::vec u(n * (n - 1) / 2);
        arma::vec &b = start_point;
        arma::vec old_b;
        // arma::mat old_z(n,n);
        // arma::mat old_u(n,n);

        z.zeros();
        // old_z.zeros();
        u.zeros();
        // old_u.zeros();

        int cnt = 0;
        do
        {
            old_b = b;
            cnt++;
            arma::vec z_row_sums(n);
            arma::vec z_col_sums(n);
            arma::vec u_row_sums(n);
            arma::vec u_col_sums(n);
            tri_sums(u, u_col_sums, u_row_sums, n);
            tri_sums(z, z_col_sums, z_row_sums, n);
            // beta subproblem: O(n)
            // TODO: vectorize
            for (int i = 0; i < n; i++)
            {
                double part1 = z_row_sums(i) + u_row_sums(i);
                double part2 = z_col_sums(i) + u_col_sums(i);
                b(i)         = ((y(i) + n * y_bar) + part1 - part2) / (n + 1);
            }
            // u and z subproblems: O(n(n-1)/2)
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    // z

                    int k              = tri_idx(i, j, n);
                    double to_be_thres = b(i) - b(j) - u(k);
                    double scale       = (1 - l * w(k) / std::abs(to_be_thres));  // TODO: vectorize
                    if (scale < 0)
                    {
                        scale = 0;
                    };
                    z(k) = scale * to_be_thres;
                    // u
                    u(k) = u(k) + (z(k) - (b(i) - b(j)));
                }
            }
            if (acc)
            {
                MoMALogger::error("Not provided yet");
                // double alpha = (1 + std::sqrt(old_alpha)) / 2;
                // z += (old_alpha / alpha) * (z - old_z);
                // old_z = z;
                // u += (old_alpha / alpha) * (u - old_u);
                // old_u = u;
                // // tri_momentum(z,old_z,old_alpha / alpha,n);
                // // tri_momentum(u,old_u,old_alpha / alpha,n);
                // old_alpha = alpha;
            }
        } while (arma::norm(old_b - b, 2) / arma::norm(old_b, 2) > prox_eps && cnt < MAX_IT);

        // Check convergence
        if (cnt == MAX_IT)
        {
            MoMALogger::warning("No convergence in unordered fusion lasso prox (ADMM).");
        }
        else
        {
            MoMALogger::debug("ADMM converges after iter: ") << cnt;
        }
        // TODO: shrink stepsize, as is done in the paper
        return b;
    }
    else
    {
        // AMA
        // Reference: Algorithm 3 Fast AMA in
        // Splitting Methods for Convex Clustering
        // Eric C. Chi and Kenneth Lange†
        int n = x.n_elem;
        // Choosing nu is not trivial. See Proposition 4.1 in
        // Splitting Methods for Convex Clustering, Chi and Range

        MoMALogger::debug("Running AMA.");
        // Find the degrees of nodes
        arma::vec deg(n);
        deg.zeros();
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (w(tri_idx(i, j, n)) > 0)
                {
                    deg(i)++;
                    deg(j)++;
                }
            }
        }

        // Find the degree of edges, which are the sums of their nodes.
        int max_edge_deg = -1;
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (w(tri_idx(i, j, n)) > 0)
                {
                    max_edge_deg = std::max(double(deg(i) + deg(j)), (double)max_edge_deg);
                }
            }
        }
        double nu    = 1.0 / (std::min(n, max_edge_deg));
        arma::vec &u = start_point;
        arma::vec lambda(n * (n - 1) / 2);
        // Initialze

        lambda.zeros();
        double old_alpha = 1;
        arma::vec old_lambda(n * (n - 1) / 2);
        old_lambda.zeros();
        arma::vec old_u;
        int cnt = 0;

        // Start iterating
        do
        {
            cnt++;
            // lambda subproblem
            for (int i = 0; i < n; i++)
            {
                for (int j = i + 1; j < n; j++)
                {
                    int k     = tri_idx(i, j, n);
                    lambda(k) = lambda(k) - nu * (u(i) - u(j));
                    // project onto the interval [ -w_ij*lambda_ij, w_ij*lambda_ij ]
                    if (std::abs(lambda(k)) > l * w(k))
                    {
                        lambda(k) = l * w(k) * lambda(k) / std::abs(lambda(k));
                    }
                }
            }

            // u subproblem
            old_u = u;
            arma::vec lambda_row_sums(n);
            arma::vec lambda_col_sums(n);
            tri_sums(lambda, lambda_col_sums, lambda_row_sums, n);
            for (int i = 0; i < n; i++)
            {
                double part1 = lambda_row_sums(i);
                double part2 = lambda_col_sums(i);
                u(i)         = x(i) + part1 - part2;
            }
            if (acc)
            {  // Momemtum step
                double alpha = (1 + std::sqrt(old_alpha)) / 2;
                tri_momentum(lambda, old_lambda, (old_alpha - 1.0) / alpha, n);
                old_alpha = alpha;
            }
        } while (arma::norm(u - old_u, 2) / arma::norm(old_u, 2) > prox_eps && cnt < MAX_IT);

        // Check convergence
        if (cnt == MAX_IT)
        {
            MoMALogger::warning("No convergence in unordered fusion lasso prox (AMA).");
        }
        else
        {
            MoMALogger::debug("AMA converges after iter: ") << cnt;
        }
        return u;
    }
}

// A handle class
ProxOp::ProxOp(Rcpp::List prox_arg_list, int dim)
{
    const std::string &s  = Rcpp::as<std::string>(prox_arg_list["P"]);
    double gamma          = Rcpp::as<double>(prox_arg_list["gamma"]);
    const arma::vec group = Rcpp::as<arma::vec>(prox_arg_list["group"]);
    double lambda2        = Rcpp::as<double>(prox_arg_list["lambda2"]);
    const arma::mat w     = Rcpp::as<arma::mat>(prox_arg_list["w"]);
    bool ADMM             = Rcpp::as<bool>(prox_arg_list["ADMM"]);
    bool acc              = Rcpp::as<bool>(prox_arg_list["acc"]);
    double prox_eps       = Rcpp::as<double>(prox_arg_list["prox_eps"]);
    int l1tf_k            = Rcpp::as<int>(prox_arg_list["l1tf_k"]);
    bool nonneg           = Rcpp::as<bool>(prox_arg_list["nonneg"]);
    if (s.compare("NONE") == 0)
    {
        p = new NullProx();
    }
    else if (s.compare("LASSO") == 0)
    {
        if (nonneg)
        {
            p = new NonNegativeLasso();
        }
        else
        {
            p = new Lasso();
        }
    }
    else if (s.compare("SCAD") == 0)
    {
        if (nonneg)
        {
            p = new NonNegativeSCAD(gamma);
        }
        else
        {
            p = new SCAD(gamma);
        }
    }
    else if (s.compare("MCP") == 0)
    {
        if (nonneg)
        {
            p = new NonNegativeMCP(gamma);
        }
        else
        {
            p = new MCP(gamma);
        }
    }
    else if (s.compare("SLOPE") == 0)
    {
        if (nonneg)
        {
            MoMALogger::error("Non-negative SLOPE is not implemented!");
        }
        else
        {
            p = new SLOPE(dim);
        }
    }
    else if (s.compare("GRPLASSO") == 0)
    {
        if (group.n_elem != dim)
        {
            MoMALogger::error("Wrong dimension: length(group) != dim(x).");
        }
        if (nonneg)
        {
            p = new NonNegativeGrpLasso(group);
        }
        else
        {
            p = new GrpLasso(group);
        }
    }
    else if (s.compare("ORDEREDFUSED") == 0)
    {
        if (nonneg)
        {
            MoMALogger::error("Non-negative ordered fused lasso is not implemented!");
        }
        else
        {
            p = new OrderedFusedLasso();
        }
    }
    else if (s.compare("ORDEREDFUSEDDP") == 0)
    {
        if (nonneg)
        {
            MoMALogger::error("Non-negative ordered fused lasso is not implemented!");
        }
        else
        {
            p = new OrderedFusedLassoDP();
        }
    }
    else if (s.compare("SPARSEFUSEDLASSO") == 0)
    {
        if (nonneg)
        {
            MoMALogger::error("Non-negative sparse fused lasso is not implemented!");
        }
        else
        {
            p = new SparseFusedLasso(lambda2);
        }
    }
    else if (s.compare("L1TRENDFILTERING") == 0)
    {
        if (nonneg)
        {
            MoMALogger::error("Non-negative L1 linear trend filtering is not implemented!");
        }
        else
        {
            p = new L1TrendFiltering(dim, l1tf_k);  // now support any order of difference matrix
        }
    }
    else if (s.compare("UNORDEREDFUSION") == 0)
    {
        if (w.n_rows != dim || w.n_cols != dim)
        {
            MoMALogger::error("Wrong dimension: dim(weight matrix) != dim(x).");
        }
        if (nonneg)
        {
            MoMALogger::error("Non-negative unordered fusion lasso is not implemented!");
        }
        else
        {
            p = new Fusion(w, ADMM, acc, prox_eps);
        }
    }
    else
    {
        MoMALogger::error(
            "Your sparse penalty is not provided by us/specified by you! Use "
            "`NullProx` by default: ")
            << s;
    }
}

arma::vec ProxOp::operator()(const arma::vec &x, double l)
{
    return (*p)(x, l);
}

int ProxOp::df(const arma::vec &x)
{
    return (*p).df(x);
}
