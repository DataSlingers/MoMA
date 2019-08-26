SFLDA <- R6::R6Class("SFLDA",
    private = list(
        check_input_index = TRUE,
        private_get_mat_by_index = function(alpha_x = 1, alpha_y = 1, lambda_x = 1, lambda_y = 1) {
            # private functions can be called only by
            # internal functions
            private$check_input_index <- FALSE
            res <- self$get_mat_by_index(
                alpha_x = alpha_x,
                alpha_y = alpha_y,
                lambda_x = lambda_x,
                lambda_y = lambda_y
            )
            private$check_input_index <- TRUE
            return(res)
        },
        private_X_project = function(newX, ...,
                                             alpha_x = 1, alpha_y = 1, lambda_x = 1, lambda_y = 1, rank = 1) {
            private$check_input_index <- FALSE
            res <- self$X_project(
                newX = newX,
                alpha_x = alpha_x,
                alpha_y = alpha_y,
                lambda_x = lambda_x,
                lambda_y = lambda_y,
                rank = rank
            )
            private$check_input_index <- TRUE
            return(res)
        },
        private_error_if_extra_arg = function(..., is_missing) {
            is_fixed <- self$fixed_list
            if (any(is_fixed == TRUE & is_missing == FALSE)) {
                param_str_list <- c("alpha_x", "alpha_y", "lambda_x", "lambda_y")
                output_para <- is_missing == FALSE & is_fixed == TRUE

                moma_error(
                    paste0(
                        "Invalid index: ",
                        paste(param_str_list[output_para], collapse = ", "),
                        ". Do not specify indexes of parameters ",
                        "i) that are chosen by BIC, or ",
                        "ii) that are not specified during initialization of the SFCCA object, or ",
                        "iii) that are scalars during initialization of the SFCCA object."
                    )
                )
            }
        },
        private_error_if_not_indices = function(...,
                                                        alpha_x, alpha_y, lambda_x, lambda_y) {
            error_if_not_finite_numeric_scalar(alpha_x)
            error_if_not_finite_numeric_scalar(alpha_y)
            error_if_not_finite_numeric_scalar(lambda_x)
            error_if_not_finite_numeric_scalar(lambda_y)

            error_if_not_wholenumber(alpha_x)
            error_if_not_wholenumber(alpha_y)
            error_if_not_wholenumber(lambda_x)
            error_if_not_wholenumber(lambda_y)
        }
    ),
    public = list(
        center_X = NULL,
        scale_X = NULL,
        grid_result = NULL,
        Omega_x = NULL,
        Omega_y = NULL,
        x_sparsity = NULL,
        y_sparsity = NULL,
        rank = NULL,
        alpha_x = NULL,
        alpha_y = NULL,
        lambda_x = NULL,
        lambda_y = NULL,
        select_scheme_list = NULL,
        pg_settings = NULL,
        n = NULL,
        px = NULL,
        py = NULL, # LDA_SPECIAL_PART, number of groups
        X = NULL,
        Y_factor = NULL, # LDA_SPECIAL_PART
        x_coln = NULL,
        x_rown = NULL,
        y_coln = NULL, # LDA_SPECIAL_PART
        fixed_list = NULL,
        initialize = function(X, ..., Y_factor, # LDA_SPECIAL_PART
                                      center = TRUE, scale = FALSE,
                                      x_sparsity = empty(), y_sparsity = empty(), lambda_x = 0, lambda_y = 0, # lambda_x/_y is a vector or scalar
                                      Omega_x = NULL, Omega_y = NULL, alpha_x = 0, alpha_y = 0, # so is alpha_x/_y
                                      pg_settings = moma_pg_settings(),
                                      select_scheme_list = list(
                                          select_scheme_alpha_x = SELECTION_SCHEME[["grid"]],
                                          select_scheme_alpha_y = SELECTION_SCHEME[["grid"]],
                                          select_scheme_lambda_x = SELECTION_SCHEME[["grid"]],
                                          select_scheme_lambda_y = SELECTION_SCHEME[["grid"]]
                                      ),
                                      max_bic_iter = 5,
                                      rank = 1,
                                      deflation_scheme = DEFLATION_SCHEME["LDA"]) {
            chkDots(...)
            # Step 1: check ALL arguments
            # Step 1.1: lambdas and alphas
            error_if_not_valid_parameters(alpha_x)
            error_if_not_valid_parameters(alpha_y)
            error_if_not_valid_parameters(lambda_x)
            error_if_not_valid_parameters(lambda_y)

            self$alpha_x <- alpha_x
            self$alpha_y <- alpha_y
            self$lambda_x <- lambda_x
            self$lambda_y <- lambda_y

            # Step 1.2: matrix
            # LDA_SPECIAL_PART
            if (!is_factor(Y_factor)) {
                moma_error("`Y` must be a factor")
            }
            # `model.matrix(~y-1)` turns a factor y into
            # an indicator matrix. `Y` is NOT stored in the
            # R6 object. Instead `Y_factor` is kept.
            X <- as.matrix(X)
            Y <- model.matrix(~ Y_factor - 1) / sqrt(length(Y_factor)[1]) # LDA_SPECIAL_PART
            error_if_not_valid_data_matrix(X)
            if (dim(X)[1] != dim(Y)[1]) {
                moma_error("`X` and `Y_factor` must have the same number of samples.")
            }

            n <- dim(X)[1]
            px <- dim(X)[2]
            py <- dim(Y)[2] # number of groups   # LDA_SPECIAL_PART
            X <- scale(X, center = center, scale = scale)

            cen_X <- attr(X, "scaled:center")
            sc_X <- attr(X, "scaled:scale")
            if (any(sc_X == 0)) {
                moma_error("cannot rescale a constant/zero column to unit variance")
            }

            self$center_X <- cen_X %||% FALSE
            self$scale_X <- sc_X %||% FALSE
            self$n <- n
            self$py <- py
            self$px <- px
            self$X <- X
            self$Y_factor <- Y_factor # LDA_SPECIAL_PART
            self$x_coln <- colnames(X) %||% paste0("Xcol_", seq_len(px))
            self$x_rown <- rownames(X) %||% paste0("Xrow_", seq_len(n))
            self$y_coln <- levels(Y_factor) # LDA_SPECIAL_PART

            # Step 1.3: sparsity
            error_if_not_of_class(x_sparsity, "_moma_sparsity_type")
            error_if_not_of_class(y_sparsity, "_moma_sparsity_type")
            self$x_sparsity <- x_sparsity
            self$y_sparsity <- y_sparsity

            # Step 1.4: PG loop settings
            error_if_not_of_class(pg_settings, "moma_pg_settings")
            self$pg_settings <- pg_settings

            # Step 1.5: smoothness
            Omega_x <- check_omega(Omega_x, alpha_x, px)
            Omega_y <- check_omega(Omega_y, alpha_y, py)
            self$Omega_x <- Omega_x
            self$Omega_y <- Omega_y

            # Step 1.6: check selection scheme string
            # `select_scheme_list` will be passed to C++ functions
            # Note `select_scheme_list` here follows _u/_v naming convention
            error_if_not_valid_select_scheme_list(select_scheme_list, uv_naming = FALSE)

            parameter_length_list <- vapply(FUN = length, list(
                self$alpha_x,
                self$alpha_y,
                self$lambda_x,
                self$lambda_y
            ), integer(1))

            self$select_scheme_list <- select_scheme_list
            self$fixed_list <- get_fixed_indicator_list(select_scheme_list, parameter_length_list, uv_naming = FALSE)

            # Step 1.7: check rank
            # TODO: check that `rank` < min(rank(X), rank(Y))
            # w.r.t to certain numeric precision
            if (!inherits(rank, "numeric") ||
                !is.wholenumber(rank) ||
                rank <= 0 ||
                rank > min(px, py, n)) { # LDA_SPECIAL_PART
                moma_error("`rank` should be a positive integer smaller than the minimum-dimension of the data matrix.")
            }
            self$rank <- rank

            # Step 2: pack all arguments in a list
            # WARNING
            # _x/_y convention on R side but _u/_v on C++ side
            algo_settings_list <- c(
                list(
                    X = X,
                    Y = Y, # LDA_SPECIAL_PART
                    lambda_u = lambda_x,
                    lambda_v = lambda_y,
                    # smoothness
                    alpha_u = alpha_x,
                    alpha_v = alpha_y,
                    rank = rank
                ),
                list(
                    Omega_u = Omega_x,
                    Omega_v = Omega_y,
                    prox_arg_list_u = add_default_prox_args(x_sparsity),
                    prox_arg_list_v = add_default_prox_args(y_sparsity)
                ),
                pg_settings,
                list(
                    select_scheme_alpha_u = select_scheme_list$select_scheme_alpha_x,
                    select_scheme_alpha_v = select_scheme_list$select_scheme_alpha_y,
                    select_scheme_lambda_u = select_scheme_list$select_scheme_lambda_x,
                    select_scheme_lambda_v = select_scheme_list$select_scheme_lambda_y
                ),
                list(
                    max_bic_iter = max_bic_iter
                ),
                list(
                    deflation_scheme = deflation_scheme # LDA_SPECIAL_PART
                )
            )
            # make sure we explicitly specify ALL arguments
            if (length(setdiff(
                names(algo_settings_list),
                names(formals(cca))
            )) != 0) {
                moma_error("Incomplete arguments in SFLDA::initialize.")
            }

            # Step 3: call the fucntion
            self$grid_result <- do.call(
                cca,
                algo_settings_list
            )
        },

        get_mat_by_index = function(..., alpha_x = 1, alpha_y = 1, lambda_x = 1, lambda_y = 1) {
            chkDots(...)

            # they should be of length 1
            private$private_error_if_not_indices(
                alpha_x = alpha_x,
                alpha_y = alpha_y,
                lambda_x = lambda_x,
                lambda_y = lambda_y
            )
            # A "fixed" parameter should not be specified
            # at all (this is a bit stringent, can be improved later).
            # "Fixed" parameters are those
            # i) that are chosen by BIC, or
            # ii) that are not specified during initialization of the SFLDA object, or
            # iii) that are scalars as opposed to vectors during initialization of the SFLDA object.

            # When `get_mat_by_index` is called internally
            # we skip the input checking
            if (private$check_input_index) {
                is_missing <- list(missing(alpha_x), missing(alpha_y), missing(lambda_x), missing(lambda_y))
                private$private_error_if_extra_arg(is_missing = is_missing)
            }

            n <- self$n
            px <- self$px
            py <- self$py
            rank <- self$rank

            U <- matrix(0, nrow = px, ncol = rank)
            V <- matrix(0, nrow = py, ncol = rank)
            d <- vector(mode = "numeric", length = rank)

            chosen_lambda_x <- vector(mode = "numeric", length = rank)
            chosen_alpha_x <- vector(mode = "numeric", length = rank)
            chosen_lambda_y <- vector(mode = "numeric", length = rank)
            chosen_alpha_y <- vector(mode = "numeric", length = rank)

            for (i in (1:self$rank)) {
                rank_i_result <- get_5Dlist_elem(self$grid_result,
                    alpha_u_i = alpha_x,
                    lambda_u_i = lambda_x,
                    alpha_v_i = alpha_y,
                    lambda_v_i = lambda_y, rank_i = i
                )[[1]]

                U[, i] <- rank_i_result$u$vector
                V[, i] <- rank_i_result$v$vector
                d[i] <- rank_i_result$d

                chosen_lambda_x[i] <- rank_i_result$u$lambda
                chosen_alpha_x[i] <- rank_i_result$u$alpha
                chosen_lambda_y[i] <- rank_i_result$v$lambda
                chosen_alpha_y[i] <- rank_i_result$v$alpha
            }


            dimnames(V) <-
                list(self$y_coln, paste0("PC", seq_len(rank)))
            dimnames(U) <-
                list(self$x_coln, paste0("PC", seq_len(rank)))
            return(list(
                X_PC_loadings = U,
                Y_group_scores = V,
                d = d,
                chosen_lambda_x = chosen_lambda_x,
                chosen_lambda_y = chosen_lambda_y,
                chosen_alpha_x = chosen_alpha_x,
                chosen_alpha_y = chosen_alpha_y
            ))
        },

        print = function() {
            selection_list_str <- lapply(self$select_scheme_list, function(x) {
                if (x == 0) {
                    return("grid search")
                }
                else if (x == 1) {
                    return("BIC search")
                }
            })

            cat("An <SFLDA> object containing solutions to the following settings\n")
            cat("Rank: ", self$rank, "\n")
            cat("Penalty and selection:\n")

            cat(paste0("alpha_x: ", selection_list_str[1], ", range: "))
            cat(self$alpha_x, "\n")
            cat(paste0("alpha_y: ", selection_list_str[2], ", range: "))
            cat(self$alpha_y, "\n")
            cat(paste0("lambda_x: ", selection_list_str[3], ", range: "))
            cat(self$lambda_x, "\n")
            cat(paste0("lambda_y: ", selection_list_str[4], ", range: "))
            cat(self$lambda_y, "\n")
        },

        X_project = function(newX, ...,
                                     alpha_x = 1, alpha_y = 1, lambda_x = 1, lambda_y = 1, rank = 1) {
            chkDots(...)
            # check indexes
            if (private$check_input_index) {
                is_missing <- list(missing(alpha_x), missing(alpha_y), missing(lambda_x), missing(lambda_y))
                private$private_error_if_extra_arg(is_missing = is_missing)
            }

            if (rank > self$rank) {
                moma_error("Invalid `rank` in SFLDA::left_project.")
            }

            private$private_error_if_not_indices(
                alpha_x = alpha_x,
                alpha_y = alpha_y,
                lambda_x = lambda_x,
                lambda_y = lambda_y
            )

            X_PC_loadings_rank_k <- private$private_get_mat_by_index(
                alpha_x = alpha_x,
                alpha_y = alpha_y,
                lambda_x = lambda_x,
                lambda_y = lambda_y
            )$X_PC_loadings[, 1:rank]


            # newX should be uncencter and unscaled.
            # check new X has same colnames
            if (length(dim(newX)) != 2L) {
                moma_error("'newX' must be a matrix or data frame")
            }

            if (dim(newX)[2] != self$px) {
                moma_error(
                    paste0(
                        "`newX` is incompatible with orignal data. ",
                        "It must have ", self$px, " columns."
                    )
                )
            }


            scaled_data <- scale(newX, self$center_X, self$scale_X)
            result <- project(scaled_data, X_PC_loadings_rank_k)
            colnames(result) <- paste0("PC", seq_len(rank))

            return(list(
                scaled_data = scaled_data,
                proj_data = result
            ))
        }
    )
)


#' The Deflation Scheme for LDA
#'
#' In \code{MoMA} one deflation scheme is provided for LDA.
#'
#' Let \eqn{X} be a data matrix (properly scaled and centered), and \eqn{Y}
#' be the indicator matrix showing which group a sample belongs to.
#' \eqn{X} and \eqn{Y} should have the same number of columns. The penalized LDA problem is formulated as
#'
#' \eqn{ \min_{u,v} \, u^T X^T Y v + \lambda_u P_u(u) + \lambda_v P_v(v)  }
#'
#' \eqn{ \text{s.t. } \| u \|_{I+\alpha_u \Omega_u} \leq 1, \| v \|_{I + \alpha_v \Omega_v} \leq 1.  }
#'
#' In the discussion below, let \eqn{u,v} be the solution to the above problem.
#' Let \eqn{c_x = Xu, c_y = Yv}. The deflation scheme is as follow:
#'
#' \eqn{X \leftarrow  { X } -  { c_x } \left(  { c_x } ^ { T }  { c_x } \right) ^ { - 1 }  { c_x } ^ { T }  { X }
#' = ( I - { c_x } \left(  { c_x } ^ { T }  { c_x } \right) ^ { - 1 }  { c_x } ^ { T } )X,}
#'
#' \eqn{ Y \text{   remains unchanged.}}.
#'
#' @references De Bie T., Cristianini N., Rosipal R. (2005) Eigenproblems
#' in Pattern Recognition. In: Handbook of Geometric Computing. Springer, Berlin, Heidelberg
#' @name LDA_deflation


NULL

#' Sparse and functional LDA
#'
#' \code{moma_sflda} creates an \code{SFLDA} R6 object and returns it.
#' @param Y_factor A factor representing which group a sample belongs to.
#' @inheritParams moma_sfcca
#' @name moma_sflda
#' @export
moma_sflda <- function(X, ..., Y_factor,
                       center = TRUE, scale = FALSE,
                       x_sparse = moma_empty(), y_sparse = moma_empty(),
                       x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                       pg_settings = moma_pg_settings(),
                       max_bic_iter = 5,
                       rank = 1) {
    chkDots(...)
    error_if_not_of_class(x_sparse, "moma_sparsity_type")
    error_if_not_of_class(y_sparse, "moma_sparsity_type")
    error_if_not_of_class(x_smooth, "moma_smoothness_type")
    error_if_not_of_class(y_smooth, "moma_smoothness_type")


    return(SFLDA$new(
        X,
        Y_factor = Y_factor,
        center = center, scale = scale,
        # sparsity
        x_sparsity = x_sparse$sparsity_type,
        y_sparsity = y_sparse$sparsity_type,
        lambda_x = x_sparse$lambda,
        lambda_y = y_sparse$lambda,
        # smoothness
        Omega_x = x_smooth$Omega,
        Omega_y = y_smooth$Omega,
        alpha_x = x_smooth$alpha,
        alpha_y = y_smooth$alpha,
        pg_settings = pg_settings,
        # Map strings to encoding
        select_scheme_list = list(
            select_scheme_alpha_x = SELECTION_SCHEME[[x_smooth$select_scheme]],
            select_scheme_alpha_y = SELECTION_SCHEME[[y_smooth$select_scheme]],
            select_scheme_lambda_x = SELECTION_SCHEME[[x_sparse$select_scheme]],
            select_scheme_lambda_y = SELECTION_SCHEME[[y_sparse$select_scheme]]
        ),
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform one-way sparse LDA
#'
#' \code{moma_slda} is a function for performing one-way sparse LDA.
#' @export
#' @describeIn moma_sflda a function for performing one-way sparse LDA
moma_slda <- function(X, ..., Y_factor,
                      center = TRUE, scale = FALSE,
                      x_sparse = moma_empty(), y_sparse = moma_empty(),
                      #    x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                      pg_settings = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1) {
    chkDots(...)
    is_x_penalized <- !missing(x_sparse)
    is_y_penalized <- !missing(y_sparse)
    if (!is_x_penalized && !is_y_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (is_x_penalized && is_y_penalized) {
        moma_error("Please use `moma_twslda` if both sides are penalized.")
    }

    return(moma_sflda(
        X = X,
        Y_factor = Y_factor,
        center = center, scale = scale,
        x_sparse = x_sparse, y_sparse = y_sparse,
        # x_smooth = x_smooth, y_smooth = y_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
    # moma_error("Not implemented: SLDA")
}


#' Perform two-way sparse LDA
#'
#' \code{moma_twslda} is a function for performing two-way sparse LDA.
#' @export
#' @describeIn moma_sflda a function for performing two-way sparse LDA
moma_twslda <- function(X, ..., Y_factor,
                        center = TRUE, scale = FALSE,
                        x_sparse = moma_empty(), y_sparse = moma_empty(),
                        #    x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                        pg_settings = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1) {
    chkDots(...)
    is_x_penalized <- !missing(x_sparse)
    is_y_penalized <- !missing(y_sparse)
    if (!is_x_penalized && !is_y_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (is_x_penalized != is_y_penalized) {
        moma_warning("Please use `moma_slda` if only one side is penalized.")
    }

    return(moma_sflda(
        X = X,
        Y_factor = Y_factor,
        center = center, scale = scale,
        x_sparse = x_sparse, y_sparse = y_sparse,
        # x_smooth = x_smooth, y_smooth = y_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform one-way functional LDA
#'
#' \code{moma_flda} is a function for performing one-way functional LDA.
#' @export
#' @describeIn moma_sflda a function for performing one-way functional LDA
moma_flda <- function(X, ..., Y_factor,
                      center = TRUE, scale = FALSE,
                      #    x_sparse = moma_empty(), y_sparse = moma_empty(),
                      x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                      pg_settings = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1) {
    chkDots(...)
    is_x_penalized <- !missing(x_smooth)
    is_y_penalized <- !missing(y_smooth)
    if (!is_x_penalized && !is_y_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (is_x_penalized && is_y_penalized) {
        moma_error("Please use `moma_twflda` if both sides are penalized.")
    }

    return(moma_sflda(
        X = X,
        Y_factor = Y_factor,
        center = center, scale = scale,
        # x_sparse = x_sparse, y_sparse = y_sparse,
        x_smooth = x_smooth, y_smooth = y_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform two-way functional LDA
#'
#' \code{moma_twflda} is a function for performing two-way functional LDA.
#' @export
#' @describeIn moma_sflda a function for performing two-way functional LDA
moma_twflda <- function(X, ..., Y_factor,
                        center = TRUE, scale = FALSE,
                        #    x_sparse = moma_empty(), y_sparse = moma_empty(),
                        x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                        pg_settings = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1) {
    chkDots(...)
    is_x_penalized <- !missing(x_smooth)
    is_y_penalized <- !missing(y_smooth)
    if (!is_x_penalized && !is_y_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (!is_x_penalized || !is_y_penalized) {
        moma_warning("Please use `moma_flda` if only one side is penalized.")
    }

    return(moma_sflda(
        X = X,
        Y_factor = Y_factor,
        center = center, scale = scale,
        # x_sparse = x_sparse, y_sparse = y_sparse,
        x_smooth = x_smooth, y_smooth = y_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}
