moma_svd <- function(
                     X,
                     u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                     Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                     pg_settings = moma_pg_settings(),
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
    # "_moma_sparsity_type" includes all penalty types
    if (!inherits(u_sparsity, "_moma_sparsity_type") || !inherits(v_sparsity, "_moma_sparsity_type")) {
        moma_error(
            "Sparse penalty should be of class ",
            sQuote("_moma_sparsity_type"),
            ". Try using, for example, `u_sparsity = lasso()`."
        )
    }

    # PG loop settings
    if (!inherits(pg_settings, "moma_pg_settings")) {
        moma_error(
            "pg_settings penalty should be of class ",
            sQuote("moma_pg_settings"),
            ". Try using, for example, `pg_settings = moma_pg_settings(MAX_ITER=1e+4)`."
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
        pg_settings
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
