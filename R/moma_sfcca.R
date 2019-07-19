SFCCA <- R6::R6Class("SFCCA",
    private = list(
        check_input_index = TRUE,
        private_get_mat_by_index = function(alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1) {
            # private functions can be called only by
            # internal functions
            private$check_input_index <- FALSE
            res <- self$get_mat_by_index(
                alpha_u, alpha_v, lambda_u, lambda_v
            )
            private$check_input_index <- TRUE
            return(res)
        }
    ),
    public = list(
        center_X = NULL,
        scale_X = NULL,
        center_Y = NULL,
        scale_Y = NULL,
        grid_result = NULL,
        Omega_u = NULL, # maybe Omega_x
        Omega_v = NULL, # maybe Omega_y
        u_sparsity = NULL, # maybe x_sparsity
        v_sparsity = NULL, # maybe y_sparsity
        rank = NULL,
        alpha_u = NULL,
        alpha_v = NULL,
        lambda_u = NULL,
        lambda_v = NULL,
        selection_scheme_list = NULL,
        pg_setting = NULL,
        n = NULL,
        px = NULL,
        py = NULL,
        X = NULL,
        Y = NULL,
        fixed_list = NULL,
        initialize = function(X, Y, ...,
                                      center = TRUE, scale = FALSE,
                                      u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                                      Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                                      pg_setting = moma_pg_settings(),
                                      selection_scheme_str = "gggg",
                                      max_bic_iter = 5,
                                      rank = 1) {
            chkDots(...) # avoid typo in argument names, e.g. y = ...
            # Step 1: check ALL arguments
            # Step 1.1: lambdas and alphas
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
            self$alpha_u <- alpha_u
            self$alpha_v <- alpha_v
            self$lambda_u <- lambda_u
            self$lambda_v <- lambda_v

            # Step 1.2: matrix
            X <- as.matrix(X)
            if (any(!is.finite(X))) {
                moma_error("X must not have NaN, NA, or Inf.")
            }
            Y <- as.matrix(Y)
            if (any(!is.finite(Y))) {
                moma_error("`Y` must not have NaN, NA, or Inf.")
            }

            if (dim(X)[1] != dim(Y)[1]) {
                moma_error("`X` and `Y` must have the same number of rows.")
            }

            n <- dim(X)[1]
            px <- dim(X)[2]
            py <- dim(Y)[2]
            X <- scale(X, center = center, scale = scale)
            Y <- scale(Y, center = center, scale = scale)

            cen_X <- attr(X, "scaled:center")
            cen_Y <- attr(Y, "scaled:center")
            sc_X <- attr(X, "scaled:scale")
            sc_Y <- attr(Y, "scaled:scale")
            if (any(sc_X == 0) || any(sc_Y == 0)) {
                moma_error("cannot rescale a constant/zero column to unit variance")
            }

            self$center_X <- cen_X %||% FALSE
            self$center_Y <- cen_Y %||% FALSE
            self$scale_X <- sc_X %||% FALSE
            self$scale_Y <- sc_Y %||% FALSE
            self$n <- dim(X)[1]
            self$py <- py
            self$px <- px
            self$X <- X
            self$Y <- Y

            # Step 1.3: sparsity
            if (!inherits(u_sparsity, "moma_sparsity") || !inherits(v_sparsity, "moma_sparsity")) {
                moma_error(
                    "Sparse penalty should be of class ",
                    sQuote("moma_sparsity"),
                    ". Try using, for example, `u_sparsity = lasso()`."
                )
            }
            self$u_sparsity <- u_sparsity
            self$v_sparsity <- v_sparsity

            # Step 1.4: PG loop settings
            if (!inherits(pg_setting, "moma_pg_settings")) {
                moma_error(
                    "pg_setting penalty should be of class ",
                    sQuote("moma_pg_settings"),
                    ". Try using, for example, `pg_setting = moma_pg_settings(MAX_ITER=1e+4)`."
                )
            }
            self$pg_setting <- pg_setting

            # Step 1.5: smoothness
            Omega_u <- check_omega(Omega_u, alpha_u, px)
            Omega_v <- check_omega(Omega_v, alpha_v, py)
            self$Omega_u <- Omega_u
            self$Omega_v <- Omega_v

            # Step 1.6: check selection scheme string
            # "g" stands for grid search, "b" stands for BIC
            if (!inherits(selection_scheme_str, "character") || nchar(selection_scheme_str) != 4 ||
                !all(strsplit(selection_scheme_str, split = "")[[1]] %in% c("b", "g"))) {
                moma_error(
                    "Invalid selection_scheme_str ", selection_scheme_str,
                    ". It should be a four-char string containing only 'b' or 'g'."
                )
            }

            # turn "b"/"g" to 1/0
            selection_scheme_list <- list(
                selection_criterion_alpha_u = 0,
                selection_criterion_alpha_v = 0,
                selection_criterion_lambda_u = 0,
                selection_criterion_lambda_v = 0
            )
            fixed_list <- list(
                # "Fixed" parameters are those
                # i) that are chosen by BIC, or
                # ii) that are not specified during initialization of the SFCPA object, or
                # iii) that are scalars as opposed to vectors during initialization of the SFCPA object.
                is_alpha_u_fixed = FALSE,
                is_alpha_v_fixed = FALSE,
                is_lambda_u_fixed = FALSE,
                is_lambda_v_fixed = FALSE
            )

            parameter_list <- list(
                self$alpha_u,
                self$alpha_v,
                self$lambda_u,
                self$lambda_v
            )
            for (i in 1:4) {
                selection_scheme_list[[i]] <- ifelse(substr(selection_scheme_str, i, i) == "g", 0, 1)
                fixed_list[[i]] <- substr(selection_scheme_str, i, i) == "b" || length(parameter_list[[i]]) == 1
            }
            self$selection_scheme_list <- selection_scheme_list
            self$fixed_list <- fixed_list

            # Step 1.7: check rank
            is.wholenumber <-
                function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
            if (!inherits(rank, "numeric") || !is.wholenumber(rank) || rank <= 0
            || rank > min(px, py)) {
                moma_error("`rank` should be a positive integer smaller than the rank of the data matrix.")
            }
            self$rank <- rank

            # Step 2: pack all arguments in a list
            algo_settings_list <- c(
                list(
                    X = X,
                    Y = Y,
                    lambda_u = lambda_u,
                    lambda_v = lambda_v,
                    # smoothness
                    alpha_u = alpha_u,
                    alpha_v = alpha_v,
                    rank = rank
                ),
                list(
                    Omega_u = Omega_u,
                    Omega_v = Omega_v,
                    prox_arg_list_u = add_default_prox_args(u_sparsity),
                    prox_arg_list_v = add_default_prox_args(v_sparsity)
                ),
                pg_setting,
                selection_scheme_list,
                list(
                    max_bic_iter = max_bic_iter
                ),
                list(
                    deflation_scheme = DEFLATION_SCHEME["CCA"]
                )
            )
            # make sure we explicitly specify ALL arguments
            if (length(setdiff(
                names(algo_settings_list),
                names(formals(cca))
            )) != 0) {
                moma_error("Incomplete arguments in SFCCA::initialize.")
            }

            # Step 3: call the fucntion
            self$grid_result <- do.call(
                cca,
                algo_settings_list
            )
        }
    )
)
