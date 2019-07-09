MOMA_EMPTYMAT <- matrix()
MOMA_EMPTYVEC <- vector(mode = "numeric")
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


# This function checks the validity of Omega and alpha
check_omega <- function(Omega, alpha, n) {

    # check if Omega is a matrix
    if (!is.matrix(Omega) && !is.null(Omega)) {
        moma_error("Omage_u/v is not a matrix.")
    }

    # store them as sparse matrices using the package Matrix
    if (length(alpha) == 1 && alpha == 0) {
        # discard the Omega matrix specified by users
        Omega <- diag(n)
    }
    else if (is.null(Omega)) {
        # The user wants smooth penalty
        # but does not specify Omega matrix
        Omega <- second_diff_mat(n)
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
    }
    return(Omega)
}

second_diff_mat <- function(n) {
    return(crossprod(diff(diag(n))))
}

moma_svd <- function(
                     X,
                     u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                     Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                     pg_setting = moma_pg_settings(),
                     k = 1, # number of pairs of singular vecters
                     select = c("gridsearch", "nestedBIC")) {
    if (!inherits(alpha_u, c("numeric", "integer")) ||
        !inherits(alpha_v, c("numeric", "integer")) ||
        !inherits(lambda_u, c("numeric", "integer")) ||
        !inherits(lambda_v, c("numeric", "integer"))) {
        moma_error(paste0(
            "All penalty levels (",
            sQuote("lambda_u"), ", ",
            sQuote("lambda_v"), ", ",
            sQuote("alpha_u"), ", ",
            sQuote("alpha_v"),
            ") must be numeric."
        ))
    }

    select <- match.arg(select)
    all_para <- c(alpha_u, alpha_v, lambda_u, lambda_v)

    # verify all alphas and lambdas are positive numbers
    if (any(all_para < 0) || any(!is.finite(all_para))) {
        moma_error(
            "All penalty levels (",
            sQuote("lambda_u"), ", ",
            sQuote("lambda_v"), ", ",
            sQuote("alpha_u"), ", ",
            sQuote("alpha_v"),
            ") must be non-negative numeric."
        )
    }

    # from scalar to vector
    alpha_u <- as.vector(alpha_u)
    alpha_v <- as.vector(alpha_v)
    lambda_u <- as.vector(lambda_u)
    lambda_v <- as.vector(lambda_v)

    if (!is.matrix(X)) {
        moma_error("X must be a matrix.")
    }
    if (any(!is.finite(X))) {
        moma_error("X must not have NaN, NA, or Inf.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]

    # If all of alpha_u, alpha_v, lambda_u, lambda_v are
    # a number, we just solve ONE MoMA problem.
    is_multiple_para <- length(alpha_u) > 1 ||
        length(alpha_v) > 1 ||
        length(lambda_u) > 1 ||
        length(lambda_v) > 1

    # k must be 1 if alpha_u/v or lambda_u/v is of vector form
    if (is_multiple_para && k != 1) {
        moma_error("We don't support a range of parameters in finding a rank-k svd")
    }

    # Sparsity arguments
    # "moma_sparsity" includes all penalty types, including fused lasso
    # group lasso and so on.
    if (!inherits(u_sparsity, "moma_sparsity") || !inherits(v_sparsity, "moma_sparsity")) {
        moma_error(
            "Sparse penalty should be of class ",
            sQuote("moma_sparsity"),
            ". Try using, for example, `u_sparsity = lasso()`."
        )
    }

    # PG loop settings
    if (!inherits(pg_setting, "moma_pg_settings")) {
        moma_error(
            "pg_setting penalty should be of class ",
            sQuote("moma_pg_settings"),
            ". Try using, for example, `pg_setting = moma_pg_settings(MAX_ITER=1e+4)`."
        )
    }

    # Pack all argument into a list
    # First we check the smoothness term argument.
    algo_settings_list <- c(
        list(
            X = X,
            lambda_u = lambda_u,
            lambda_v = lambda_v,
            # smoothness
            alpha_u = alpha_u,
            alpha_v = alpha_v,
            rank = k
        ),
        # Penalties
        list(
            Omega_u = check_omega(Omega_u, alpha_u, n),
            Omega_v = check_omega(Omega_v, alpha_v, p),
            prox_arg_list_u = add_default_prox_args(u_sparsity),
            prox_arg_list_v = add_default_prox_args(v_sparsity)
        ),
        pg_setting
    )

    if (is_multiple_para) {
        if (select == "gridsearch") {
            a <- do.call("cpp_moma_grid_search", algo_settings_list)
            class(a) <- "moma_svd_grid"
            return(a)
        }
        else if (select == "nestedBIC") {
            a <- do.call("cpp_moma_criterion_search", algo_settings_list)
            class(a) <- "moma_svd_nestedBIC"
            return(a)
        }
        else {
            moma_error("Wrong parameter selection methods!")
        }
    }
    else {
        return(do.call("cpp_moma_multi_rank", algo_settings_list))
    }
}
