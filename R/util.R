`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs)) {
        lhs
    } else {
        rhs
    }
}

is_fin_numeric_scalar <- function(x) {
    (length(x) == 1L) && is.numeric(x) && (!is.na(x)) && is.finite(x)
}
is_valid_parameters <- function(x) {
    length(x) >= 1 && all(sapply(x, is_fin_numeric_scalar))
}
is_valid_select_str <- function(x) {
    is.character(x) && nchar(x) == 1 && x %in% c("b", "g")
}
# Check whether `x` is a boolean value
is_logical_scalar <- function(x) {
    return(is.logical(x) && (length(x) == 1) && !is.na(x))
}
# Credit: clustRviz
is_square <- function(x) {
    is.matrix(x) && (NROW(x) == NCOL(x))
}

MOMA_EMPTYMAT <- matrix() # the default `w` arguemnt for unordered fusion
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

# This function checks the validity of Omega and alpha
check_omega <- function(Omega, alpha, n) {

    # check if Omega is a matrix or a NULL
    if (!is.matrix(Omega) && !is.null(Omega)) {
        moma_error("Omega_u/v is not a matrix.")
    }

    # TODO: store them as sparse matrices using the package Matrix
    if (length(alpha) == 1 && alpha == 0) {
        # discard the Omega matrix specified by users
        Omega <- diag(n) # TODO: should not overwrite
    }
    else if (is.null(Omega)) {
        # The user wants smooth penalty
        # but does not specify Omega matrix
        Omega <- second_diff_mat(n) # TODO: should not overwrite
    }
    else {
        # At this point, users have specified an Omega and
        # non-zero penalty levels explicitly

        # Check validity of Omega if users speicify both alpha and Omega
        if (dim(Omega)[1] != dim(Omega)[2]) {
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

second_diff_mat <- function(n) {
    return(crossprod(diff(diag(n))))
}
