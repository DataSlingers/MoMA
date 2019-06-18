#include "moma_prox.h"

/*
 * L1 Trend Filtering
 */

// A second difference matrix
// [0 0 ... 1 -2 1 ... 0 0]
// [[Rcpp::export]]
arma::mat l1tf_diff_mat(int m, int k)
{
    if (k < 0)
    {
        MoMALogger::error("k should be non-negative integer.");
    }
    if (m < k + 2)
    {
        MoMALogger::error("A difference matrix should have more columns.");
    }

    arma::vec coef(k + 2);  // Note for k = 0, D = [... 1, -1 ...]
                            // k = 1, D = [... 1, -2, 1 ...]
                            // k = 2, D = [... 1, -3, 3, -1, ...]
                            // The (k+1)-th row of the Pascal triangle
    coef(0)     = std::pow(-1, k);
    int flag    = -coef(0);
    coef(k + 1) = -1;
    for (int i = 1; i <= (k + 1) / 2; i++)
    {
        coef(i)         = -1 * coef(i - 1) * (k - i + 2) / i;
        coef(k + 1 - i) = flag * coef(i);
    }

    arma::mat D = arma::zeros<arma::mat>(m - 1 - k, m);
    for (int i = 0; i < m - 1 - k; i++)
    {
        for (int j = 0; j < k + 2; j++)
        {
            D(i, i + j) = coef(j);
        }
    }

    return D;
}

L1TrendFiltering::L1TrendFiltering(int n, int i_k)
{
    if (i_k >= 2)
    {
        MoMALogger::message("TF with higher-than-second difference matrix is not well tested yet.");
    }
    if (n == -1)
    {
        MoMALogger::error("Class L1TrendFiltering needs to know dimension of the problem.");
    }

    D = l1tf_diff_mat(n, i_k);

    k = i_k;

    MoMALogger::debug(
        "Initializing a L1 linear trend filtering proximal operator object of "
        "degree ")
        << k << ".";
}

L1TrendFiltering::~L1TrendFiltering()
{
    MoMALogger::debug("Releasing a L1 linear trend filtering proximal operator object");
}

// Find the sum of dual residual and centering residual
double res(const arma::mat &DDt,
           const arma::vec &Dy,
           const arma::vec &mu1,
           const arma::vec &mu2,
           const arma::vec &f1,
           const arma::vec &f2,
           double t,
           const arma::vec &nu)
{
    arma::vec res_cent = -(mu1 % f1) - (mu2 % f2) - 1 / t;
    // NOTE: this approximates formulea 16 in the paper
    arma::vec res_dual = DDt * nu - Dy + mu1 - mu2;
    return arma::as_scalar(arma::norm(res_cent) + arma::norm(res_dual));
}

double init_stepsize(const arma::vec &mu, const arma::vec &dmu, double step)
{
    if (step <= 0)
    {
        MoMALogger::error("step should > 0 in init_stepsize.");
    }

    arma::uvec idx = (dmu < 0);

    if (arma::sum(idx) > 0)
    {
        // the largest step size that does not negate any elements of mu
        // i.e., max(step) s.t. mu + step * dmu > 0 for all elements
        arma::vec prop = -(idx % (mu / dmu));

        double max_legit_ss = step;
        for (int i = 0; i < mu.n_elem; i++)
        {
            if (prop(i) > 0 && prop(i) < max_legit_ss)
            {
                max_legit_ss = 0.99 * prop(i);
            }
        }
        return max_legit_ss;
    }
    else
    {
        return step;
    }
}

arma::vec L1TrendFiltering::operator()(const arma::vec &x, double l)
{
    int n = x.n_elem;
    int m = n - 1 - k;  // # of rows of the diff mat
    if (D.n_cols != n || D.n_rows != m)
    {
        MoMALogger::error("Error in initialzing difference matrix.");
    }
    const arma::vec &y = x;

    // Commonly used mat and vec
    arma::mat DDt = D * D.t();
    arma::vec Dy  = D * y;

    // Step size
    double step;

    // Initialzation with a strictly feasible point
    // Ref: l1 Trend Filtering by Seung-Jean Kim, Kwangmoo Koh
    // Stephen Boyd and Dimitry Gorinevsky
    arma::vec nu  = arma::zeros<arma::vec>(m);  // The dual variable, y - D^t*nu is what we need
    arma::vec mu1 = arma::ones<arma::vec>(m);   // Multipliers to solve the dual problem
    arma::vec mu2 = arma::ones<arma::vec>(m);

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
    for (; iter < MAX_ITER; iter++)
    {
        // Surrogate duality gap
        // Ref: Primal-Dual Interior-Point Methods,
        // Ryan Tibshirani, Convex Optimization 10-725/36-725
        // Powerpoint presentation page 10
        // http://www.cs.cmu.edu/~pradeepr/convexopt/Lecture_Slides/primal-dual.pdf
        gap = -(sum(mu1 % f1 + mu2 % f2));
        if (gap < prox_eps)
        {
            break;
        }

        // Ref: PPT page 11
        // NOTE: this is different from the `l1tf` C implementation
        t = 4 * m / gap;

        arma::mat part1 = DDt - arma::diagmat(mu1 / f1 + mu2 / f2);  // A band matrix
        arma::vec part2 = -DDt * nu + Dy + (1 / t) / f1 - (1 / t) / f2;

        // Update directions
        arma::vec dnu =
            arma::solve(part1,
                        part2);  // RcppArmadillo optimizes with band matrix automatically
        arma::vec dmu1 = -(mu1 + ((1 / t) + dnu % mu1) / f1);
        arma::vec dmu2 = -(mu2 + ((1 / t) + dnu % mu2) / f2);

        // Choose step size by back tracking
        // Initialize with the largest stepsize not making mu's negative
        step = init_stepsize(mu1, dmu1, 1);
        step = init_stepsize(mu2, dmu2, step);

        for (int iter_bt = 0; iter_bt < MAX_BT_ITER; iter_bt++)
        {
            new_nu  = nu + step * dnu;
            new_mu1 = mu1 + step * dmu1;
            new_mu2 = mu2 + step * dmu2;

            // new inequality constraints
            new_f1 = new_nu - l;
            new_f2 = -new_nu - l;

            // check step size, exit if all of the conditions are met
            // 1. f1<0 and f2<0
            // 2. loss function should be less than a reference linear function
            if (arma::sum(new_f1 > 0) > 0 || arma::sum(new_f2 > 0) > 0 ||
                res(DDt, Dy, new_mu1, new_mu2, new_f1, new_f2, t, new_nu) >
                    (1 - alpha * step) * res(DDt, Dy, mu1, mu2, f1, f2, t, nu))
            {
                step *= beta;
            }
            else
            {
                break;
            }
        }
        MoMALogger::debug("Picking step = ") << step << ".";
        // Update variable with selected stepsize
        nu  = new_nu;
        mu1 = new_mu1;
        mu2 = new_mu2;
        f1  = new_f1;
        f2  = new_f2;
    }

    if (iter == MAX_ITER)
    {
        MoMALogger::info(
            "No convergence in L1 linear trend filtering solver. Surrogate duality "
            "(gap,iter) = ")
            << "(" << gap << "," << iter << ")"
            << ".";
    }
    arma::vec Dymm = Dy - (mu1 - mu2);
    double dres =
        0.5 * arma::as_scalar(Dymm.t() * arma::solve(DDt, Dymm)) + l * arma::sum(mu1 + mu2);
    arma::vec Dtnu = D.t() * nu;
    double pres    = -0.5 * arma::as_scalar(Dtnu.t() * Dtnu) + arma::as_scalar(Dy.t() * nu);
    MoMALogger::debug("Primal loss = ")
        << pres << " , dual loss = " << dres << " , gap = " << pres - dres
        << " , surrogate gap = " << gap << " , iter = " << iter << ".";
    return y - D.t() * nu;
}

int L1TrendFiltering::df(const arma::vec &x)
{
    if (k == 0)
    {
        MoMALogger::error("Please use fused lasso instead.");
    }
    else if (k == 1)
    {
        // df = number of knots + k + 1
        // Ref:
        // Table 2 of Tibshirani, Robert.
        // "Regression shrinkage and selection via the lasso: a retrospective."
        // Journal of the Royal Statistical Society: Series B (Statistical
        // Methodology) 73.3 (2011): 273-282. D = [... −1 2 −1 ...]

        if (x.n_elem < 3)
        {
            MoMALogger::error("dim(x) must be larger than 2.");
        }
        int knots = 0;
        // Count the number of knots (changes in slope in the case of linear trend
        // filtering)
        for (int i = 2; i < x.n_elem; i++)
        {
            if (std::abs(x(i - 2) - 2 * x(i - 1) + x(i)) > 1e-10)
            {
                printf("%d", i);
                knots++;
            }
        }
        return knots + k + 1;
    }
    else if (k == 2)
    {
        // df = number of knots + k + 1
        // Ref:
        // Table 2 of Tibshirani, Robert.
        // "Regression shrinkage and selection via the lasso: a retrospective."
        // Journal of the Royal Statistical Society: Series B (Statistical
        // Methodology) 73.3 (2011): 273-282. D = [...   1   -3    3   -1  ...]

        if (x.n_elem < 4)
        {
            MoMALogger::error("dim(x) must be larger than 3.");
        }
        int knots = 0;
        // Count the number of knots (changes in second derivative in the case of
        // quadratic trend filtering)
        for (int i = 3; i < x.n_elem; i++)
        {
            if (std::abs(x(i - 3) - 3 * x(i - 2) + 3 * x(i - 1) - x(i)) > 1e-10)
            {
                knots++;
            }
        }

        return knots + k + 1;
    }
    else
    {
        MoMALogger::error("Error in L1TrendFiltering::df: Invalid k.");
        return 0;
    }
}
