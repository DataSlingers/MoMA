#' Sparsity-inducing penalty in \code{MoMA}
#' In the package \code{MoMA}, we support the following sparsity-inducing
#' penalty functions.
#' \itemize{
#'   \item{\code{\link{lasso}}} TODO
#'   \item{\code{\link{mcp}}} TODO
#'   \item{\code{\link{scad}}} TODO
#'   \item{\code{\link{slope}}} TODO
#'   \item{\code{\link{grplasso}}} TODO
#'   \item{\code{\link{fusedlasso}}} TODO
#'   \item{\code{\link{l1tf}}} TODO
#'   \item{\code{\link{spfusedlasso}}} TODO
#'   \item{\code{\link{cluster}}} TODO
#' }
#' @name moma_sparsity
NULL

empty <- function() {
    arglist <- list()
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' LASSO
#'
#' Use this function to set the penalty function as lasso
#' \deqn{\lambda \sum \| x_{i} \| ,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param non_negative a boolean value. Set \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return a \code{moma_sparsity} object, which is a list containing \code{non_negative}
#'
#' @examples
#' lasso(non_negative = FALSE)
#' @export
lasso <- function(..., non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    if (!is_logical_scalar(non_negative)) {
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    arglist <- list(nonneg = non_negative, P = "LASSO")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' MCP (minimax concave penalty)
#'
#' Use this function to set the penalty function as MCP
#' \deqn{\lambda P (x; \gamma),}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}, \eqn{P} is
#' determined by \eqn{\gamma}. See Zhang, Cun-Hui. "Nearly unbiased variable
#' selection under minimax concave penalty." The Annals of statistics 38.2 (2010): 894-942.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param gamma non-convexity. Must be larger than 1.
#' @param non_negative a boolean value. Set to \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return a \code{moma_sparsity} object, which is a list containing \code{non_negative}
#' and \code{gamma}.
#'
#' @examples
#' mcp(gamma = 3, non_negative = FALSE)
#' @export
mcp <- function(..., gamma = 3, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    if (!is_logical_scalar(non_negative)) {
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if (gamma <= 1) {
        moma_error(
            "Non-convexity parameter of MCP (",
            sQuote("gamma"),
            ") must be larger than 1."
        )
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "MCP")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' SCAD
#'
#' Use this function to set the penalty function as SCAD
#' \deqn{\lambda P (x; \gamma) ,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}, \eqn{P} is
#' determined by \eqn{\gamma}. See Fan, Jianqing, and Runze Li. "Variable selection
#'  via nonconcave penalized likelihood and its oracle properties." Journal of
#'  the American statistical Association 96.456 (2001): 1348-1360.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param gamma non-convexity. Must be larger than 2.
#' @param non_negative a boolean value. Set to \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return a \code{moma_sparsity} object, which is a list containing \code{non_negative}
#' and \code{gamma}.
#'
#' @examples
#' scad(gamma = 3.7, non_negative = FALSE)
#' @export
scad <- function(..., gamma = 3.7, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    if (!is_logical_scalar(non_negative)) {
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if (gamma <= 2) {
        moma_error(
            "Non-convexity parameter of SCAD (",
            sQuote("gamma"),
            ") must be larger than 2."
        )
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "SCAD")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' SLOPE
#'
#' Use this function to set the penalty function as SLOPE - Sorted L-One Penalized Estimation
#' \deqn{\lambda P (x) = \lambda \sum _ { i = 1 } ^ { n } \lambda _ { i } | x | _ { ( i ) } .}
#' where \eqn{\lambda_i = \Phi ^ { - 1 } ( 1 - q _ { i } ) ,  q _ { i } = i \cdot q / 2 p, q = 0.05.}
#' Here \eqn{q} is the false discovery rate (FDR).
#' See Bogdan, Malgorzata, et al. "SLOPE - adaptive variable selection via convex optimization."
#' The annals of applied statistics 9.3 (2015): 1103.
#'
#' @return a \code{moma_sparsity} object, which contains a list containing the string "SLOPE".
#'
#' @examples
#' slope()
#' @export
slope <- function() {
    arglist <- list(P = "SLOPE")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' Group LASSO
#'
#' Use this function to set the penalty function as group lasso
#' \deqn{\lambda \sum_{g \in Group} \| x_g \|,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}, \eqn{\|x_g\|} is
#' the vector comprised of elements of \eqn{x} picked out by indeces set \eqn{g}.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param g a vector of integer or characters, or a factor itself. It gets transformed
#' to factor eventually to indicate grouping.
#' @param non_negative a boolean value. Set to \code{TRUE} to add non-negativity
#' constraint.
#'
#' @return a \code{moma_sparsity} object, which is a list containing \code{non_negative}
#' and \code{g}.
#'
#' @examples
#' # This sets every three adjacent parameters as a group.
#' grplasso(g = rep(1:10, each = 3), non_negative = FALSE)
#' @export
grplasso <- function(..., g, non_negative = FALSE) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }

    if (!is_logical_scalar(non_negative)) {
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if (!(inherits(g, c("character", "numeric", "factor", "integer")))) {
        moma_error("Please provide a vector as an indicator of grouping.")
    }
    arglist <- list(group = as.factor(g), P = "GRPLASSO", nonneg = non_negative)
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' Fused lasso
#'
#' Use this function to set the penalty function as fused lasso
#' \deqn{\lambda \sum \| x_{i} - x_{i-1} \|,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}.
#' @param ... Forces users to specify all arguments by name.
#' @param algo a string being either "path" or "dp". Defaults to "path". Partial matching
#' is supported. Two solving algorithms
#' are provided. When "path" is chosen, the algorithm by
#' Hoefling, H. (2010) is used. When "dp" is chosen, the algorithm by Johnson, N. A. (2013) is used.
#'
#' @references Hoefling, H. (2010). A path algorithm
#' for the fused lasso signal approximator. Journal of Computational and Graphical
#'  Statistics, 19(4), 984-1006, doi: 10.1198/jcgs.2010.09208.
#' @references Johnson, N. A. (2013). A dynamic programming algorithm for the
#' fused lasso and l 0-segmentation. Journal of Computational and Graphical
#' Statistics, 22(2), 246-260, doi: 10.1080/10618600.2012.681238.
#' @return a \code{moma_sparsity} object, which is an empty list.
#'
#' @examples
#' fusedlasso()
#' @export
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
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' L1 trend filtering
#'
#' Use this function to set the penalty function as l1 trend filtering. An
#' important special case is when \eqn{k=1}. In this case the penalty
#' term becomes
#' \deqn{\lambda \sum \| x_{i-1} - 2x_{i} + x_{i+1} \|,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}.
#' For other values of \eqn{k} please refer to the following table:
#' \tabular{llll}{
#' k=0                \tab k=1              \tab k=2                \tab ... \cr
#' piecewise constant \tab peicewise linear \tab piecewise quadratic \tab ...
#' }
#' The general formula of the penalty term for \eqn{k \in N} can be found in
#' Tibshirani, Ryan J. "Adaptive piecewise polynomial estimation via trend
#' filtering." The Annals of Statistics 42.1 (2014): 285-323.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param l1tf_k use (k+1)-difference matrix in trend filtering. Note \eqn{k = 0}
#'          implies piecewise constant, \eqn{k=1} implies piecewise linear, \eqn{k=2}
#'          piecewise quadratic etc.
#'
#' @return a \code{moma_sparsity} object, which is an empty list.
#'
#' @examples
#' l1tf(l1tf_k = 1)
#' @export
l1tf <- function(..., l1tf_k = 1) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    # l1 linear trend filtering
    arglist <- list(P = "L1TRENDFILTERING", l1tf_k = l1tf_k)
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' Sparse fused lasso
#'
#' Use this function to set the penalty function as sparse fused lasso
#' \deqn{\lambda_1 \sum \| x_{i} - x_{i-1} \| + \lambda_2 \sum \|x_{i} \| ,}
#' where \eqn{\lambda_} is set by \code{lambda_u/v} in the function \code{moma_svd}, and \eqn{\lambda_2}
#' is specified in here.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param lambda2 the level of penalty on the absolute values of the coefficients
#'
#' @return a \code{moma_sparsity} object, which is a list containing the value of \code{lambda_2}.
#'
#' @examples
#' spfusedlasso(lambda2 = 2)
#' @export
spfusedlasso <- function(..., lambda2) {
    if (length(list(...)) != 0) {
        moma_error("Please specify the correct argument by name.")
    }
    arglist <- list(P = "SPARSEFUSEDLASSO", lambda2 = lambda2)
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

#' Cluster penalty
#'
#' Use this function to set the penalty function as
#' \deqn{\lambda \sum w_{ij} \| x_{i} - x_{j} \|,}
#' where \eqn{\lambda} is set by \code{lambda_u/v} in the function \code{moma_svd}.
#'
#' @param ... Forces users to specify all arguments by name.
#' @param w a symmetric square matrix. \code{w[i][j]} is the \eqn{w_{ij}} described above.
#' @param ADMM a boolean value. Set to \code{TRUE} to use ADMM, set to \code{FALSE} to use AMA.
#' @param acc a boolean value. Set to \code{TRUE} to use the accelereated version of the algorithm.
#' Currently we support accelerated AMA only.
#' @param eps a small numeric value. The precision used when solving the proximal operator.
#'
#' @return a \code{moma_sparsity} object, which is a list containing the values of \code{w},
#' \code{ADMM}, \code{acc} and \code{eps}.
#'
#' @examples
#' cluster(w = matrix(rep(1, 9), 3), ADMM = FALSE, acc = FALSE, eps = 1e-10)
#' @export
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
    class(arglist) <- "moma_sparsity"
    return(arglist)
}


#' Algorithm settings for solving a penalzied SVD problem
#'
#' To find an (approximate) solution to a penalized SVD (Singular Value Decomposition) problem is to solve two
#' penalized regression problems iteratively (outer loop). Each penalized regression (inner loop)
#' is solved using one of the three algorithms: ISTA (Iterative Shrinkage-Thresholding Algorithm),
#' FISTA (Fast Iterative Shrinkage-Thresholding Algorithm) and
#' One-step ISTA (an approximated version of ISTA).
#' @param ... to force users to specify arguments by names
#' @param EPS precision for outer loop
#' @param MAX_ITER the maximum number of iterations for outer loop
#' @param EPS_inner precision for inner loop
#' @param MAX_ITER_inner the maximum number of iterations for inner loop
#' @param solver a string in \code{c("ista", "fista", "onestepista")}.
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
        if (!is_valid_parameters(lambda)) {
            moma_error(sQuote("lambda"), " is not valid: ", lambda)
        }

        # Step 3: check select_scheme
        if (!is_valid_select_str(select_scheme)) {
            moma_error(
                sQuote("select_scheme"), " is not valid: ", select_scheme,
                ". It should be either `g` or `b`."
            )
        }

        # step 4: return
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
moma_lasso <- create_moma_sparsity_func(lasso)
moma_mcp <- create_moma_sparsity_func(mcp)
moma_scad <- create_moma_sparsity_func(scad)
moma_slope <- create_moma_sparsity_func(slope)
moma_grplasso <- create_moma_sparsity_func(grplasso)
moma_fusedlasso <- create_moma_sparsity_func(fusedlasso)
moma_l1tf <- create_moma_sparsity_func(l1tf)
moma_spfusedlasso <- create_moma_sparsity_func(spfusedlasso)
moma_cluster <- create_moma_sparsity_func(cluster)

# What this function does now is just wrap three arguement
# into a list.
# TODO: `Omega` could be a user-defined function
moma_smoothness <- function(Omega = NULL, ..., alpha = 0, select_scheme = "g") {

    # Step 2: check lambda
    if (!is_valid_parameters(alpha)) {
        moma_error(sQuote("alpha"), " is not valid: ", alpha)
    }

    # Step 3: check select_scheme
    if (!is_valid_select_str(select_scheme)) {
        moma_error(
            sQuote("select_scheme"), " is not valid: ", select_scheme,
            ". It should be either `g` or `b`."
        )
    }

    a <- list(
        Omega = Omega,
        alpha = alpha,
        select_scheme = select_scheme
    )
    class(a) <- "moma_smoothness_type"
    return(a)
}
