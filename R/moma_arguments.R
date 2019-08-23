#' Sparsity-inducing penalty in \code{MoMA}
#'
#' In the package \code{MoMA}, we support the following sparsity-inducing
#' penalty functions.
#' \itemize{
#'   \item{\code{\link{moma_lasso}}}: sparsity
#'   \item{\code{\link{moma_mcp}}}: non-convex sparsity
#'   \item{\code{\link{moma_scad}}}: non-convex sparsity
#'   \item{\code{\link{moma_slope}}}: sparsity
#'   \item{\code{\link{moma_grplasso}}}: group-wise sparsity
#'   \item{\code{\link{moma_fusedlasso}}}: piecewise constant, or ordered fusion
#'   \item{\code{\link{moma_spfusedlasso}}}: sparsity and piece-wise constant
#'   \item{\code{\link{moma_l1tf}}}: piecewise polynomial (default to piecewise linear)
#'   \item{\code{\link{moma_cluster}}}: unordered fusion
#' }
#' These functions specify the value of the \code{u_sparse,v_sparse} arguments in the
#'  \code{moma_*pca} series of functions, and the \code{x_sparse,y_sparse} arguments
#' in the \code{moma_*cca} and \code{moma_*lda} series of functions.
#'
#' All functions
#' above share two common parameters: \code{lambda} and \code{select_scheme}, which are
#' described in the Arguments section.
#' @name moma_sparsity_options
#' @param ... Force users to specify arguments by names.
#' @param lambda A vector containing penalty values
#' @param select_scheme A char being either "b" (nested BIC search) or "g" (grid search).
#'
#' MoMA provides a flexible framework for regularized multivariate analysis
#' with several tuning parameters for different forms of regularization.
#' To assist the user in selecting these parameters (\code{alpha_u},
#' \code{alpha_v}, \code{lambda_u}, \code{lambda_v}), we provide
#' two selection modes: grid search ("g") and nested BIC search ("b").
#' Grid search means we solve the problem
#' for all combinations of parameter values provided by the user.
#'
#' To explain nested BIC search, we need to look into how the algorithm runs.
#' To find an (approximate) solution to a penalized SVD (Singular Value Decomposition) problem is to solve two
#' penalized regression problems iteratively. Let's call them problem u and problem v, which give
#' improving estimates of the right singular vector, \emph{u}, and the left singular vector, \emph{v}, respectively.
#' For each regression problem, we can select the optimal parameters
#' based on BIC.
#'
#' The nested BIC search is essentially two 2-D searches. We start from SVD solutions, and then find the optimal
#' parameters for problem u, given current estimate of \emph{v}. Using the result from previous step, update
#' current estimate of \emph{u}, and then do the same thing for problem v,
#' that is, to find the optimal parameters for problem v given current estimate of \emph{u}. Repeat
#' the above until convergence or the maximal number of iterations has been reached.
#'
#' Users are welcome to refer to section 3.1: Selection of Regularization Parameters
#' in the paper cited below.
#'
#' @references G. I. Allen and M. Weylandt, "Sparse and Functional Principal
#' Components Analysis," 2019 IEEE Data Science Workshop (DSW),
#' Minneapolis, MN, USA, 2019, pp. 11-16. \doi{10.1109/DSW.2019.8755778}.
NULL

#' Introduction to selection schemes in MoMA
#'
#' Please see the description of the argument \code{select_scheme} in
#' \code{\link{moma_sparsity_options}}. The \code{select_scheme} argument
#' presents in functions listed in \code{\link{moma_sparsity_options}}, and the function
#' \code{\link{moma_smoothness}}.
#' @name select_scheme
NULL

empty <- function() {
    arglist <- list()
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

lasso <- function(..., non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    error_if_not_logical_scalar(non_negative)

    arglist <- list(nonneg = non_negative, P = "LASSO")
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

mcp <- function(..., gamma = 3, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    error_if_not_logical_scalar(non_negative)

    if (gamma <= 1) {
        moma_error(
            "Non-convexity parameter of MCP (",
            sQuote("gamma"),
            ") must be larger than 1."
        )
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "MCP")
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

scad <- function(..., gamma = 3.7, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    error_if_not_logical_scalar(non_negative)

    if (gamma <= 2) {
        moma_error(
            "Non-convexity parameter of SCAD (",
            sQuote("gamma"),
            ") must be larger than 2."
        )
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "SCAD")
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

slope <- function() {
    arglist <- list(P = "SLOPE")
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

grplasso <- function(..., g, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    error_if_not_logical_scalar(non_negative)

    if (!(inherits(g, c("character", "numeric", "factor", "integer")))) {
        moma_error("Please provide a vector as an indicator of grouping.")
    }
    arglist <- list(group = as.factor(g), P = "GRPLASSO", nonneg = non_negative)
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

fusedlasso <- function(..., algo = c("path", "dp")) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    # fused lasso

    # Two options for solving the proximal operator
    # of fused lasso are supported
    # "dp": dynamic programming
    # "path": solution path-based algorithm
    algo <- match.arg(algo)
    prox_name <- ifelse(algo == "path", "ORDEREDFUSED", "ORDEREDFUSEDDP")
    arglist <- list(P = prox_name)
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

l1tf <- function(..., l1tf_k = 1) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    # l1 linear trend filtering
    arglist <- list(P = "L1TRENDFILTERING", l1tf_k = l1tf_k)
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

spfusedlasso <- function(..., lambda2) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    arglist <- list(P = "SPARSEFUSEDLASSO", lambda2 = lambda2)
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}

cluster <- function(..., w = NULL, ADMM = FALSE,
                    acc = FALSE,
                    eps = 1e-10) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    # fused lasso
    if (!is.matrix(w) || is.null(w) || dim(w)[1] != dim(w)[2]) {
        moma_error("`w` should be a square matrix.")
    }
    if (!isSymmetric(w)) {
        moma_warning("`w` is not symmetric. Only upper triangular half is used.")
    }
    arglist <- list(
        w = w, ADMM = ADMM,
        acc = acc, prox_eps = eps,
        P = "UNORDEREDFUSION"
    )
    class(arglist) <- "_moma_sparsity_type"
    return(arglist)
}


#' Algorithm settings for solving the penalized SVD problem
#'
#' This function is used to specify the \code{pg_settings} argument
#' in the \code{moma_*pca}, \code{moma_*cca}, and \code{moma_*lda}
#' family of functions.
#'
#' To find an (approximate) solution to a penalized SVD (Singular Value Decomposition) problem is to solve two
#' penalized regression problems iteratively (outer loop). Each penalized regression (inner loop)
#' is solved using one of the three algorithms: ISTA (Iterative Shrinkage-Thresholding Algorithm),
#' FISTA (Fast Iterative Shrinkage-Thresholding Algorithm) and
#' One-step ISTA (an approximated version of ISTA).
#' @param ... To force users to specify arguments by names.
#' @param EPS Precision for outer loop.
#' @param MAX_ITER The maximum number of iterations for outer loop.
#' @param EPS_inner Precision for inner loop.
#' @param MAX_ITER_inner The maximum number of iterations for inner loop.
#' @param solver A string in \code{c("ista", "fista", "onestepista")}, representing ISTA (Iterative Shrinkage-Thresholding Algorithm),
#'              FISTA (Fast
#'              Iterative Shrinkage-Thresholding Algorithm) and One-step ISTA (an approximated
#'              version of ISTA) respectively.
#' @return A \code{moma_pg_settings} object, which is a list containing the above parameters.
#' @export
moma_pg_settings <- function(..., EPS = 1e-10, MAX_ITER = 1000,
                             EPS_inner = 1e-10, MAX_ITER_inner = 1e+5,
                             solver = c("ista", "fista", "onestepista")) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    solver <- match.arg(solver)
    arglist <- list(
        EPS = EPS, MAX_ITER = MAX_ITER,
        EPS_inner = EPS_inner, MAX_ITER_inner = MAX_ITER_inner,
        solver = toupper(solver)
    )
    class(arglist) <- "moma_pg_settings"
    return(arglist)
}

create_moma_sparsity_func <- function(f) {
    # Given f, we want to generate a new function, which
    # 1. contains all arguments in f;
    # 2. has two extra arguments `lambda` and `select_scheme`;
    # 3. returns a list that contains ( f(...), lambda = ..., select_scheme = ... ).
    aug_f <- function(..., lambda = 0, select_scheme = "g") {
        chkDots(...)

        # Step 2: check lambda
        error_if_not_valid_parameters(lambda)

        # Step 3: check select_scheme
        error_if_not_valid_select_str(select_scheme)

        select_scheme <- match_selection_scheme(select_scheme)

        # step 4: return
        # WARNING: do not define local variables before
        # `mget(ls())`
        arg_for_f <- mget(ls())
        arg_for_f <- arg_for_f[
            names(arg_for_f) %in% c("lambda", "select_scheme") == FALSE
        ] # fetch arguments for the function `f`
        a <- list(
            sparsity_type = do.call(f, arg_for_f),
            lambda = lambda,
            select_scheme = select_scheme
        )
        class(a) <- "moma_sparsity_type"
        return(a)
    }
    formals(aug_f) <- c(formals(f), formals(aug_f))
    aug_f
}

moma_empty <- create_moma_sparsity_func(empty)

#' LASSO (least absolute shrinkage and selection operator)
#'
#' Use this function to set the penalty function to LASSO
#' \deqn{\lambda \sum | x_{i} | = \lambda \| x \|_1 ,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param non_negative A Boolean value. Set \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the value of \code{non_negative}
#'
#' @references Tibshirani, Robert. "Regression Shrinkage and Selection via the Lasso."
#' Journal of the Royal Statistical Society: Series B (Methodological) 58.1 (1996): 267-288. \doi{10.1111/j.2517-6161.1996.tb02080.x}.
#' @name lasso
#' @inheritParams moma_sparsity_options
#' @export

moma_lasso <- create_moma_sparsity_func(lasso)

#' MCP (minimax concave penalty)
#'
#' Use this function to set the penalty function to MCP
#' \deqn{ P (x; \lambda, \gamma) =
#' \left\{\begin{array}{ll}{\lambda|x|-\frac{x^{2}}{2 \gamma},} & {
#' \text { if }|x| \leq \gamma \lambda} \\ {\frac{1}{2} \gamma
#' \lambda^{2},} & {\text { if }|x|>\gamma \lambda}\end{array}\right.,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#' @param ... Forces users to specify all arguments by name.
#' @param gamma Non-convexity. Must be larger than 1.
#' @param non_negative A Boolean value. Set to \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the value of \code{non_negative}
#' and \code{gamma}.
#' @references Zhang, Cun-Hui. "Nearly Unbiased Variable
#' Selection under Minimax Concave Penalty." The Annals of Statistics 38.2 (2010): 894-942. \doi{10.1214/09-AOS729}.
#'
#' @name mcp
#' @inheritParams moma_sparsity_options
#' @export
moma_mcp <- create_moma_sparsity_func(mcp)

#' SCAD (smoothly clipped absolute deviation)
#'
#' Use this function to set the penalty function to SCAD
#' \deqn{ P (x; \lambda, \gamma) = \left\{\begin{array}{ll}{
#' \lambda|x|} & {\text { if }|x| \leq \lambda} \\ {\frac{2 \gamma \lambda|x|-x^{2}-
#' \lambda^{2}}{2(\gamma-1)}} & {\text { if } \lambda<|x|<\gamma \lambda} \\
#' {\frac{\lambda^{2}(\gamma+1)}{2}} & {\text { if }|x| \geq \gamma \lambda}\end{array}\right.,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param gamma Non-convexity. Must be larger than 2.
#' @param non_negative A Boolean value. Set to \code{TRUE} to add non-negativity
#'                  constraint.
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the values of \code{non_negative}
#' and \code{gamma}.
#' @references Fan, Jianqing, and Runze Li. "Variable Selection
#'  via Nonconcave Penalized Likelihood and Its Oracle Properties." Journal of
#'  the American Statistical Association 96.456 (2001): 1348-1360. \doi{10.1198/016214501753382273}.
#' @name scad
#' @inheritParams moma_sparsity_options
#' @export
moma_scad <- create_moma_sparsity_func(scad)

#' SLOPE (sorted \eqn{\ell}-one penalized estimation)
#'
#' Use this function to set the penalty function to SLOPE - Sorted L-One Penalized Estimation
#' \deqn{\lambda P (x) = \lambda \sum _ { i = 1 } ^ { n } \lambda _ { i } | x | _ { ( i ) } .}
#' where \eqn{\lambda_i = \Phi ^ { - 1 } ( 1 - q _ { i } ) ,  q _ { i } = i \cdot q / 2 p, q = 0.05.}
#' Here \eqn{q} is the false discovery rate (FDR).
#' @references Bogdan, Malgorzata, et al. "SLOPE - Adaptive Variable Selection via Convex Optimization."
#' The Annals of Applied Statistics 9.3 (2015): 1103. \doi{10.1214/15-AOAS842}.
#'
#' @return A \code{moma_sparsity_type} object, which contains a list containing the string "SLOPE".
#' @name slope
#' @inheritParams moma_sparsity_options
#' @export
moma_slope <- create_moma_sparsity_func(slope)

#' Group LASSO
#'
#' Use this function to set the penalty function to group lasso
#' \deqn{\lambda \sum_{g \in Group} \| x_g \|_1,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below, \eqn{x_g} is
#' the vector comprised of elements of \eqn{x} picked out by the indices set \eqn{g}.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param g A vector of integer or characters, or a factor itself. It gets transformed
#' to factor eventually to indicate grouping.
#' @param non_negative A Boolean value. Set it to \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the values of \code{non_negative}
#' and \code{g}.
#' @name group_lasso
#' @references Yuan, Ming, and Yi Lin. "Model Selection and Estimation in Regression
#' with Grouped Variables." Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology) 68.1 (2006): 49-67. \doi{10.1111/j.1467-9868.2005.00532.x}.
#' @inheritParams moma_sparsity_options
#' @export
moma_grplasso <- create_moma_sparsity_func(grplasso)

#' Fused LASSO
#'
#' Use this function to set the penalty function to fused lasso
#' \deqn{\lambda \sum | x_{i} - x_{i-1} |,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#' @param ... Forces users to specify all arguments by name.
#' @param algo A string being either "path" or "dp". Defaults to "path". Partial matching
#' is supported. Two solving algorithms
#' are provided. When "path" is chosen, the algorithm by
#' Hoefling, H. (2010) is used. When "dp" is chosen, the algorithm by Johnson, N. A. (2013) is used.
#'
#' @references Tibshirani, Robert, et al. "Sparsity and Smoothness via the Fused Lasso."
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67.1 (2005): 91-108.
#' \doi{10.1111/j.1467-9868.2005.00490.x}.
#' @references Hoefling, H. (2010). A path algorithm
#' for the fused lasso signal approximator. Journal of Computational and Graphical
#'  Statistics, 19(4), 984-1006. \doi{10.1198/jcgs.2010.09208}.
#' @references Johnson, N. A. (2013). A dynamic programming algorithm for the
#' fused lasso and l 0-segmentation. Journal of Computational and Graphical
#' Statistics, 22(2), 246-260. \doi{10.1080/10618600.2012.681238}.
#' @return A \code{moma_sparsity_type} object, which is an empty list.
#'
#' @name fused_lasso
#' @inheritParams moma_sparsity_options
#' @export
moma_fusedlasso <- create_moma_sparsity_func(fusedlasso)

#' L1 trend filtering
#'
#' Use this function to set the penalty function to l1 trend filtering. An
#' important special case is when \eqn{k=1}. In this case the penalty
#' term becomes
#' \deqn{\lambda \sum | x_{i-1} - 2x_{i} + x_{i+1} |,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#'
#' The general formula of the penalty term for \eqn{k \in N} can be found in
#' the paper cited in Reference. For other values of \eqn{k} please refer to the following table:
#' \tabular{llll}{
#'   \tab \eqn{k = 0} \tab \eqn{k = 1} \tab \eqn{k = 2}
#' \cr
#'  Type of sparsity \tab piecewise constant \tab peicewise linear \tab piecewise quadratic
#' }
#'
#' @param ... Forces users to specify all arguments by name.
#' @param l1tf_k Use (\eqn{k+1})-difference matrix in trend filtering. Note \eqn{k = 0}
#'          implies piecewise constant, \eqn{k=1} implies piecewise linear, \eqn{k=2}
#'          piecewise quadratic etc.
#' @references Tibshirani, Ryan J. "Adaptive Piecewise Polynomial Estimation via Trend
#' Filtering." The Annals of Statistics 42.1 (2014): 285-323. \doi{10.1214/13-AOS1189}.
#' @references Aaditya Ramdas & Ryan J. Tibshirani (2016) Fast and Flexible ADMM
#' Algorithms for Trend Filtering, Journal of Computational and Graphical Statistics,
#' 25:3, 839-858. \doi{10.1080/10618600.2015.1054033}.
#'
#' @return A \code{moma_sparsity_type} object, which is an empty list.
#' @name l1_trend_filtering
#' @inheritParams moma_sparsity_options
#' @export
moma_l1tf <- create_moma_sparsity_func(l1tf)

#' Sparse fused LASSO
#'
#' Use this function to set the penalty function to sparse fused lasso
#' \deqn{\lambda \sum | x_{i} - x_{i-1} | + \lambda_2 \sum |x_{i} | ,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below, and \eqn{\lambda_2}
#' is specified in by the \code{lambda_2} argument.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param lambda2 A scalar. The level of penalty on the absolute values of the coefficients.
#'                Note that it remains fixed when searching over \code{lambda}, rather than
#'                changes with \code{lambda} in a way that the \code{lambda} / \code{lambda_2}
#'                ratio remains fixed (which is the defualt behavior in the package
#'                \code{glmnet}).
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the value of \code{lambda_2}.
#' @references Tibshirani, Robert, et al. "Sparsity and Smoothness via the Fused Lasso."
#' Journal of the Royal Statistical Society: Series B (Statistical Methodology) 67.1 (2005): 91-108.
#' \doi{10.1111/j.1467-9868.2005.00490.x}.
#' @name sparse_fused_lasso
#' @inheritParams moma_sparsity_options
#' @export
moma_spfusedlasso <- create_moma_sparsity_func(spfusedlasso)

#' Cluster penalty
#'
#' Use this function to set the penalty function to
#' \deqn{\lambda \sum w_{ij} | x_{i} - x_{j} |,}
#' where \eqn{\lambda} is set by the \code{lambda} argument below.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param w A symmetric square matrix. \code{w[i, j]} is the \eqn{w_{ij}} described above.
#' @param ADMM A Boolean value. Set to \code{TRUE} to use ADMM, set to \code{FALSE} to use AMA. Defaults to FALSE.
#' @param acc A Boolean value. Set to \code{TRUE} to use the accelerated version of the algorithm.
#'          Currently we support accelerated AMA only.
#' @param eps A small numeric value. The precision used when solving the proximal operator.
#'
#' @return A \code{moma_sparsity_type} object, which is a list containing the values of \code{w},
#' \code{ADMM}, \code{acc} and \code{eps}.
#'
#' @name cluster
#' @inheritParams moma_sparsity_options

#' @references Chi, Eric C., and Kenneth Lange. "Splitting Methods for Convex Clustering."
#' Journal of Computational and Graphical Statistics 24.4 (2015): 994-1013. \doi{10.1080/10618600.2014.948181}.
#' @export
moma_cluster <- create_moma_sparsity_func(cluster)


#' Smoothness-inducing Term
#'
#' This function specifies the value of the \code{u_smooth,v_smooth} arguments in the
#'  \code{moma_*pca} series of functions, and the \code{x_smooth,y_smooth} arguments
#' in the \code{moma_*cca} and \code{moma_*lda} series of functions.
#' @param Omega A matrix of appropriate size. A common choice is the second difference matrix.
#' See \code{\link{second_diff_mat}}.
#' @param alpha A vector containing penalty values
#' @inheritParams moma_sparsity_options
#' @export
moma_smoothness <- function(Omega = NULL, ..., alpha = 0, select_scheme = "g") {

    # What this function does now is just wrap three arguement
    # into a list.
    # TODO: `Omega` could be a user-defined function
    # Step 2: check lambda
    error_if_not_valid_parameters(alpha)

    # Step 3: check select_scheme
    error_if_not_valid_select_str(select_scheme)

    select_scheme <- match_selection_scheme(select_scheme)

    a <- list(
        Omega = Omega,
        alpha = alpha,
        select_scheme = select_scheme
    )
    class(a) <- "moma_smoothness_type"
    return(a)
}
