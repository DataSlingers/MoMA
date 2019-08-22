`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs)) {
        lhs
    } else {
        rhs
    }
}

is_valid_data_matrix <- function(x) {
    return(is.double(x) && all(is.finite(x)))
}
error_if_not_valid_data_matrix <- function(x) {
    nm <- deparse(substitute(x))
    if (!is_valid_data_matrix(x)) {
        moma_error(
            sQuote(nm),
            " must contain numbers only and must not have NaN, NA, or Inf."
        )
    }
}

is_finite_numeric_scalar <- function(x) {
    (length(x) == 1L) && is.numeric(x) && (!is.na(x)) && is.finite(x)
}
error_if_not_finite_numeric_scalar <- function(x) {
    nm <- deparse(substitute(x))
    if (!is_finite_numeric_scalar(x)) {
        moma_error(
            sQuote(nm),
            " must be a finite scalar."
        )
    }
}

is_valid_parameters <- function(x) {
    length(x) >= 1 && all(sapply(x, is_finite_numeric_scalar))
}
error_if_not_valid_parameters <- function(x) {
    nm <- deparse(substitute(x))
    if (!is_valid_parameters(x)) {
        moma_error(
            sQuote(nm),
            " is not a valid grid."
        )
    }
}

is_valid_select_str <- function(x) {
    is.character(x) && nchar(x) == 1 && x %in% c("b", "g")
}
error_if_not_valid_select_str <- function(x) {
    nm <- deparse(substitute(x))
    if (!is_valid_select_str(x)) {
        moma_error(
            sQuote(nm),
            " should be either `g` or `b`."
        )
    }
}
# Check whether `x` is a boolean value
is_logical_scalar <- function(x) {
    return(is.logical(x) && (length(x) == 1) && !is.na(x))
}
error_if_not_logical_scalar <- function(x) {
    nm <- deparse(substitute(x))
    if (!is_logical_scalar(x)) {
        moma_error(
            sQuote(nm),
            " should be a Boolean value."
        )
    }
}

# Credit: clustRviz
is_square <- function(x) {
    is.matrix(x) && (NROW(x) == NCOL(x))
}
is_factor <- function(a) {
    return(is.finite(a) && is.factor(a))
}

error_if_not_of_class <- function(x, cl) {
    nm <- deparse(substitute(x))
    if (!inherits(x, cl)) {
        moma_error(
            sQuote(nm),
            " must be of class ",
            sQuote(cl),
            "."
        )
    }
}

error_if_not_fourchar_bg_string <- function(x) {
    nm <- deparse(substitute(x))

    if (!inherits(x, "character") ||
        nchar(x) != 4 ||
        !all(strsplit(x, split = "")[[1]] %in% c("b", "g"))) {
        moma_error(
            sQuote(nm),
            " should be a four-char string containing only 'b's and 'g's."
        )
    }
}

MOMA_EMPTYMAT <- matrix() # the default `w` argument for unordered fusion
MOMA_EMPTYVEC <- vector(mode = "numeric") # the defautl `group` arguement for group lasso
# MOMA_DEFAULT_PROX must be consistent with
# the initializer of the C++ class ProxOp
# in `moma_prox.h`
MOMA_DEFAULT_PROX <- list(
    P = "NONE",
    gamma = 3,
    # non-negativity
    nonneg = FALSE,
    # grouping
    group = MOMA_EMPTYVEC,
    lambda2 = 0,
    # unordered fusion
    w = MOMA_EMPTYMAT,
    ADMM = FALSE,
    acc = FALSE,
    prox_eps = 1e-10,
    # trend filtering
    l1tf_k = 1
)
add_default_prox_args <- function(sparsity_type) {
    # sparsity_type: prox arguments for u and v

    # To call a C function we have to specify
    # all arguments. However, some arguments
    # are specific for a particular prox. So
    # we first assign a default arg list to
    # `df_prox_arg_list_u/_v` and
    # then update them.
    return(modifyList(MOMA_DEFAULT_PROX, sparsity_type))
}

is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
error_if_not_wholenumber <- function(x) {
    nm <- deparse(substitute(x))
    if (!is.wholenumber(x)) {
        moma_error(
            sQuote(nm),
            " must be a whole number."
        )
    }
}

# This function checks the validity of Omega and alpha
check_omega <- function(Omega, alpha, n) {

    # check if Omega is a matrix or a NULL
    if (!is.matrix(Omega) && !is.null(Omega)) {
        moma_error("Omega_u/v is not a matrix.")
    }

    ## LOGIC:
    # if alpha = 0: overwrite Omega_u to identity matrix whatever it was
    # if alpha is a grid or a non-zero scalar:
    #       if Omega missing: set to second-difference matrix
    #       else check validity

    # TODO: store them as sparse matrices using the package Matrix
    if (length(alpha) == 1 && alpha == 0) {
        # discard the Omega matrix specified by users
        Omega <- diag(n) # TODO: should not overwrite
    }
    else if (is.null(Omega)) {
        # The user wants smooth penalty
        # but does not specify Omega matrix
        Omega <- second_diff_mat(n)
        moma_info("Use the second difference matrix as the default smoothing matrix.")
    }
    else {
        # At this point, users have specified an Omega and
        # non-zero penalty levels explicitly

        # Check validity of Omega if users speicify both alpha and Omega
        if (!is_square(Omega)) {
            moma_error(
                "Omega shoud be a square matrix: nrows = ", dim(Omega)[1],
                ", ncols = ", dim(Omega)[2]
            )
        }
        if (dim(Omega)[1] != n) {
            moma_error(
                "Omega shoud be a compatible matrix. It should be of ",
                n, "x", n,
                ", but is actually ", dim(Omega)[1], "x", dim(Omega)[1]
            )
        }
        # TODO: check definiteness and symmetry of Omega
    }
    return(Omega)
}


#' Second difference matrix
#'
#' This function returns a second difference matrix of size \eqn{n}.
#'
#' @param n An integer. The size of the returned matrix.
#' @name second_diff_mat
#' @export
second_diff_mat <- function(n) {
    return(crossprod(diff(diag(n))))
}

DEFLATION_SCHEME <- c(
    PCA_Hotelling = 1,
    CCA = 2,
    LDA = 3,
    PLS = 4,
    PCA_Schur_complement = 5,
    PCA_Projection = 6
)

# project rows of X to the column
# space of V
project <- function(X, V) {
    PV <- solve(crossprod(V), t(V)) # project onto the span of V
    result <- X %*% t(PV)
    return(result)
}
