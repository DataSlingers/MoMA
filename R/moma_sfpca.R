SFPCA <- R6::R6Class("SFPCA",
    private = list(
        check_input_index = TRUE,
        private_get_mat_by_index = function(alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1) {
            # private functions can be called only by
            # internal functions
            private$check_input_index <- FALSE
            res <- self$get_mat_by_index(
                alpha_u = alpha_u,
                alpha_v = alpha_v,
                lambda_u = lambda_u,
                lambda_v = lambda_v
            )
            private$check_input_index <- TRUE
            return(res)
        },
        private_left_project = function(newX, ...,
                                                alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1, rank = 1) {
            private$check_input_index <- FALSE
            res <- self$left_project(
                newX = newX,
                alpha_u = alpha_u,
                alpha_v = alpha_v,
                lambda_u = lambda_u,
                lambda_v = lambda_v,
                rank = rank
            )
            private$check_input_index <- TRUE
            return(res)
        },
        private_error_if_not_indices = function(...,
                                                        alpha_u, alpha_v, lambda_u, lambda_v) {
            error_if_not_finite_numeric_scalar(alpha_u)
            error_if_not_finite_numeric_scalar(alpha_v)
            error_if_not_finite_numeric_scalar(lambda_u)
            error_if_not_finite_numeric_scalar(lambda_v)

            error_if_not_wholenumber(alpha_u)
            error_if_not_wholenumber(alpha_v)
            error_if_not_wholenumber(lambda_u)
            error_if_not_wholenumber(lambda_v)
        },
        private_error_if_extra_arg = function(..., is_missing) {
            is_fixed <- self$fixed_list == TRUE
            param_str_list <- c("alpha_u", "alpha_v", "lambda_u", "lambda_v")
            if (any(is_missing == FALSE & is_fixed == TRUE)) { # elementwise and
                output_para <- is_missing == FALSE & is_fixed == TRUE
                moma_error(
                    paste0(
                        "Invalid index: ",
                        paste(param_str_list[output_para], collapse = ", "),
                        ". Do not specify indexes of parameters ",
                        "i) that are chosen by BIC, or ",
                        "ii) that are not specified during initialization of the SFPCA object, or ",
                        "iii) that are scalars during initialization of the SFPCA object."
                    )
                )
            }
        },
        private_error_if_miss_arg = function(..., is_missing) {
            is_fixed <- self$fixed_list == TRUE
            param_str_list <- c("alpha_u", "alpha_v", "lambda_u", "lambda_v")
            if (any(is_missing == TRUE & is_fixed == FALSE)) {
                output_para <- is_missing == TRUE & is_fixed == FALSE
                moma_error(
                    paste0(
                        "Please spesify the following argument(s): ",
                        paste(param_str_list[output_para], collapse = ", "),
                        "."
                    )
                )
            }
        }
    ),
    public = list(
        center = NULL,
        scale = NULL,
        grid_result = NULL,
        Omega_u = NULL,
        Omega_v = NULL,
        u_sparsity = NULL,
        v_sparsity = NULL,
        rank = NULL,
        alpha_u = NULL,
        alpha_v = NULL,
        lambda_u = NULL,
        lambda_v = NULL,
        select_scheme_list = NULL,
        pg_settings = NULL,
        n = NULL,
        p = NULL,
        X = NULL,
        fixed_list = NULL,
        X_coln = NULL,
        X_rown = NULL,
        initialize = function(X, ...,
                                      center = TRUE, scale = FALSE,
                                      u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                                      Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                                      pg_settings = moma_pg_settings(),
                                      select_scheme_list = list(
                                          select_scheme_alpha_u = SELECTION_SCHEME[["grid"]],
                                          select_scheme_alpha_v = SELECTION_SCHEME[["grid"]],
                                          select_scheme_lambda_u = SELECTION_SCHEME[["grid"]],
                                          select_scheme_lambda_v = SELECTION_SCHEME[["grid"]]
                                      ),
                                      max_bic_iter = 5,
                                      rank = 1,
                                      deflation_scheme = DEFLATION_SCHEME[["PCA_Hotelling"]]) {
            chkDots(...)

            # Step 1: check ALL arguments
            # Step 1.1: lambdas and alphas
            error_if_not_valid_parameters(alpha_u)
            error_if_not_valid_parameters(alpha_v)
            error_if_not_valid_parameters(lambda_u)
            error_if_not_valid_parameters(lambda_v)

            self$alpha_u <- alpha_u
            self$alpha_v <- alpha_v
            self$lambda_u <- lambda_u
            self$lambda_v <- lambda_v

            # Step 1.2: matrix
            X <- as.matrix(X)
            error_if_not_valid_data_matrix(X)
            n <- dim(X)[1]
            p <- dim(X)[2]
            X <- scale(X, center = center, scale = scale)
            self$X_coln <- colnames(X) %||% paste0("Xcol_", seq_len(p))
            self$X_rown <- rownames(X) %||% paste0("Xrow_", seq_len(n))

            cen <- attr(X, "scaled:center")
            sc <- attr(X, "scaled:scale")
            if (any(sc == 0)) {
                moma_error("cannot rescale a constant/zero column to unit variance")
            }

            self$center <- cen %||% FALSE
            self$scale <- sc %||% FALSE
            self$n <- n
            self$p <- p
            self$X <- X
            self$X_coln <- colnames(X) %||% paste0("Xcol_", seq_len(p))
            self$X_rown <- rownames(X) %||% paste0("Xrow_", seq_len(n))

            # Step 1.3: sparsity
            error_if_not_of_class(u_sparsity, "_moma_sparsity_type")
            error_if_not_of_class(v_sparsity, "_moma_sparsity_type")
            self$u_sparsity <- u_sparsity
            self$v_sparsity <- v_sparsity

            # Step 1.4: PG loop settings
            error_if_not_of_class(pg_settings, "moma_pg_settings")
            self$pg_settings <- pg_settings

            # Step 1.5: smoothness
            # if alpha = 0: overwrite Omega_u to identity matrix whatever it was
            # if alpha is a grid or a non-zero scalar:
            #       if Omega missing: set to second-difference matrix
            #       else check validity
            Omega_u <- check_omega(Omega_u, alpha_u, n)
            Omega_v <- check_omega(Omega_v, alpha_v, p)
            self$Omega_u <- Omega_u
            self$Omega_v <- Omega_v

            # Step 1.6: check selection scheme string
            # `select_scheme_list` will be passed to C++ functions
            error_if_not_valid_select_scheme_list(select_scheme_list)

            parameter_length_list <- vapply(FUN = length, list(
                self$alpha_u,
                self$alpha_v,
                self$lambda_u,
                self$lambda_v
            ), integer(1)) # order is important

            self$select_scheme_list <- select_scheme_list
            self$fixed_list <- get_fixed_indicator_list(select_scheme_list, parameter_length_list)

            # Step 1.7: check rank
            # TODO: check that `rank` < min(rank(X), rank(Y))
            if (!inherits(rank, "numeric") ||
                !is.wholenumber(rank) ||
                rank <= 0 ||
                rank > min(p, n)) {
                moma_error(
                    sQuote("rank"),
                    " should be a positive integer smaller than the minimum-dimension of the data matrix."
                )
            }
            self$rank <- rank

            # Step 2: pack all arguments in a list
            algo_settings_list <- c(
                list(
                    X = X,
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
                pg_settings,
                select_scheme_list,
                list(
                    max_bic_iter = max_bic_iter
                ),
                list(
                    deflation_scheme = deflation_scheme
                )
            )
            # make sure we explicitly specify ALL arguments
            if (length(setdiff(
                names(algo_settings_list),
                names(formals(cpp_multirank_BIC_grid_search))
            )) != 0) {
                moma_error("Incomplete arguments in SFPCA::initialize.")
            }

            # Step 3: call the fucntion
            self$grid_result <- do.call(
                cpp_multirank_BIC_grid_search,
                algo_settings_list
            )
        },

        get_mat_by_index = function(..., alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1) {
            chkDots(...)

            private$private_error_if_not_indices(
                alpha_u = alpha_u,
                alpha_v = alpha_v,
                lambda_u = lambda_u,
                lambda_v = lambda_v
            )

            # A "fixed" parameter should not be specified
            # at all (this is a bit stringent, can be improved later).
            # "Fixed" parameters are those
            # i) that are chosen by BIC, or
            # ii) that are not specified during initialization of the SFPCA object, or
            # iii) that are scalars as opposed to vectors during initialization of the SFPCA object.

            # When `get_mat_by_index` is called internally
            # we skip the input checking
            if (private$check_input_index) {
                is_missing <- list(missing(alpha_u), missing(alpha_v), missing(lambda_u), missing(lambda_v))
                private$private_error_if_extra_arg(is_missing = is_missing)
            }

            n <- self$n
            p <- self$p
            rank <- self$rank

            U <- matrix(0, nrow = n, ncol = rank)
            V <- matrix(0, nrow = p, ncol = rank)
            d <- vector(mode = "numeric", length = rank)

            chosen_lambda_u <- vector(mode = "numeric", length = rank)
            chosen_alpha_u <- vector(mode = "numeric", length = rank)
            chosen_lambda_v <- vector(mode = "numeric", length = rank)
            chosen_alpha_v <- vector(mode = "numeric", length = rank)

            for (i in (1:self$rank)) {
                rank_i_result <- get_5Dlist_elem(self$grid_result,
                    alpha_u_i = alpha_u,
                    lambda_u_i = lambda_u,
                    alpha_v_i = alpha_v,
                    lambda_v_i = lambda_v, rank_i = i
                )[[1]]

                U[, i] <- rank_i_result$u$vector
                V[, i] <- rank_i_result$v$vector
                d[i] <- rank_i_result$d

                chosen_lambda_u[i] <- rank_i_result$u$lambda
                chosen_alpha_u[i] <- rank_i_result$u$alpha
                chosen_lambda_v[i] <- rank_i_result$v$lambda
                chosen_alpha_v[i] <- rank_i_result$v$alpha
            }


            dimnames(V) <-
                list(self$X_coln, paste0("PC", seq_len(rank)))
            dimnames(U) <-
                list(self$X_rown, paste0("PC", seq_len(rank)))
            return(list(
                U = U, V = V, d = d,
                chosen_lambda_u = chosen_lambda_u,
                chosen_lambda_v = chosen_lambda_v,
                chosen_alpha_u = chosen_alpha_u,
                chosen_alpha_v = chosen_alpha_v
            ))
        },

        interpolate = function(..., alpha_u = 0, alpha_v = 0, lambda_u = 0, lambda_v = 0, exact = FALSE) {
            chkDots(...)

            # If BIC scheme has been used for any parameters, exit.
            if (any(self$select_scheme_list != 0)) {
                moma_error("R6 object SFPCA do not support interpolation when BIC selection scheme has been used.")
            }

            # Reject inputs like alpha_u = "1" or "alpha_u" = c(1,2,3)
            error_if_not_finite_numeric_scalar(alpha_u)
            error_if_not_finite_numeric_scalar(alpha_v)
            error_if_not_finite_numeric_scalar(lambda_u)
            error_if_not_finite_numeric_scalar(lambda_v)

            # Parameters that are specified explictly is not "fixed".
            # Parameters that are "fixed" must not be specified.
            is_missing <- list(missing(alpha_u), missing(alpha_v), missing(lambda_u), missing(lambda_v))
            private$private_error_if_extra_arg(is_missing = is_missing)
            private$private_error_if_miss_arg(is_missing = is_missing)

            if (exact) {
                alpha_u <- ifelse(self$fixed_list$is_alpha_u_fixed, self$alpha_u, alpha_u)
                alpha_v <- ifelse(self$fixed_list$is_alpha_v_fixed, self$alpha_v, alpha_v)
                lambda_u <- ifelse(self$fixed_list$is_lambda_u_fixed, self$lambda_u, lambda_u)
                lambda_v <- ifelse(self$fixed_list$is_lambda_v_fixed, self$lambda_v, lambda_v)

                a <- moma_svd(
                    X = self$X,
                    # smoothness
                    u_sparsity = self$u_sparsity, v_sparsity = self$v_sparsity,
                    lambda_u = lambda_u, lambda_v = lambda_v,
                    # sparsity
                    Omega_u = self$Omega_u, Omega_v = self$Omega_v,
                    alpha_u = alpha_u, alpha_v = alpha_v,
                    pg_settings = self$pg_settings,
                    k = self$rank
                )
                return(list(U = a$u, V = a$v))
            }

            # Function `findInterval` requires sorted breakpoints
            if (is.unsorted(self$alpha_u) ||
                is.unsorted(self$alpha_v) ||
                is.unsorted(self$lambda_u) ||
                is.unsorted(self$lambda_v)) {
                moma_error("Penalty levels not sorted!")
            }

            # a bool: want to interpolate on u side?
            inter_u <- !missing(alpha_u) && !missing(lambda_u) &&
                !self$fixed_list$is_alpha_u_fixed &&
                !self$fixed_list$is_lambda_u_fixed &&
                self$fixed_list$is_alpha_v_fixed &&
                self$fixed_list$is_lambda_v_fixed

            inter_v <- !missing(alpha_v) && !missing(lambda_v) &&
                self$fixed_list$is_alpha_u_fixed &&
                self$fixed_list$is_lambda_u_fixed &&
                !self$fixed_list$is_alpha_v_fixed &&
                !self$fixed_list$is_lambda_v_fixed

            # If both of them are ture or false at the same time, exit.
            if (inter_u == inter_v) {
                moma_error("SFPCA::interpolate only supports one-sided interpolation.")
            }

            if (inter_v) {

                # test that it is in the known range
                if (alpha_v >= max(self$alpha_v) || alpha_v <= min(self$alpha_v)) {
                    moma_error("Invalid range: alpha_v.")
                }

                # find the cloest alpha
                alpha_v_i <- which.min(abs(self$alpha_v - alpha_v))

                # find the bin where lambda lies in
                if (lambda_v >= max(self$lambda_v) || lambda_v <= min(self$lambda_v)) {
                    moma_error("Invalid range: lambda_v.")
                }
                lambda_v_i_lo <- findInterval(lambda_v, self$lambda_v)
                lambda_v_i_hi <- lambda_v_i_lo + 1
                if (lambda_v_i_hi > length(self$lambda_v) || lambda_v_i_lo <= 0) {
                    moma_error("SFPCA::interpolate, error in findInterval")
                }

                result_lo <- private$private_get_mat_by_index(
                    alpha_u = 1,
                    alpha_v = alpha_v_i,
                    lambda_u = 1,
                    lambda_v = lambda_v_i_lo
                )
                result_hi <- private$private_get_mat_by_index(
                    alpha_u = 1,
                    alpha_v = alpha_v_i,
                    lambda_u = 1,
                    lambda_v = lambda_v_i_hi
                )

                U <- 0.5 * (result_lo$U + result_hi$U)
                V <- 0.5 * (result_lo$V + result_hi$V)

                # newSv <- diag(self$p) + alpha_v * self$Omega_v
                # newSu <- diag(self$n) + alpha_u * self$Omega_u

                # LvT <- chol(newSv)
                # LuT <- chol(newSu) # L^T L = newS, L is a lower triagular matrix

                return(list(U = U, V = V))
            }
            else if (inter_u) {
                if (alpha_u >= max(self$alpha_u) || alpha_u <= min(self$alpha_u)) {
                    moma_error("Invalid range: alpha_u.")
                }
                # find the cloest alpha
                alpha_u_i <- which.min(abs(self$alpha_u - alpha_u))


                # find the bin where lambda lies in
                if (lambda_u >= max(self$lambda_u) || lambda_u <= min(self$lambda_u)) {
                    moma_error("Invalid range: lambda_u.")
                }
                lambda_u_i_lo <- findInterval(lambda_u, self$lambda_u)
                lambda_u_i_hi <- lambda_u_i_lo + 1
                if (lambda_u_i_hi > length(self$lambda_u) || lambda_u_i_lo <= 0) {
                    moma_error("SFPCA::interpolate, error in findInterval.")
                }

                result_lo <- private$private_get_mat_by_index(
                    alpha_v = 1,
                    alpha_u = alpha_u_i,
                    lambda_v = 1,
                    lambda_u = lambda_u_i_lo
                )
                result_hi <- private$private_get_mat_by_index(
                    alpha_v = 1,
                    alpha_u = alpha_u_i,
                    lambda_v = 1,
                    lambda_u = lambda_u_i_hi
                )

                U <- 0.5 * (result_lo$U + result_hi$U)
                V <- 0.5 * (result_lo$V + result_hi$V)
                return(list(U = U, V = V))
            }
            else {
                moma_error("UNKNOWN.")
            }
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

            cat("An <SFPCA> object containing solutions to the following settings\n")
            cat("rank =", self$rank, "\n")
            cat("Penalty and selection:\n")

            cat(paste0("alpha_u: ", selection_list_str[1]), "\n")
            print(self$alpha_u)
            cat(paste0("alpha_v: ", selection_list_str[2]), "\n")
            print(self$alpha_v)
            cat(paste0("lambda_u: ", selection_list_str[3]), "\n")
            print(self$lambda_u)
            cat(paste0("lambda_v: ", selection_list_str[4]), "\n")
            print(self$lambda_v)
        },

        left_project = function(newX, ...,
                                        alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1, rank = 1) {

            # check indexes
            if (private$check_input_index) {
                # The order is important!
                is_missing <- list(missing(alpha_u), missing(alpha_v), missing(lambda_u), missing(lambda_v))
                private$private_error_if_extra_arg(is_missing = is_missing)
            }

            if (rank > self$rank) {
                moma_error("Invalid `rank` in SFPCA::left_project.")
            }

            private$private_error_if_not_indices(
                alpha_u = alpha_u,
                alpha_v = alpha_v,
                lambda_u = lambda_u,
                lambda_v = lambda_v
            )

            V <- private$private_get_mat_by_index(
                alpha_u = alpha_u,
                alpha_v = alpha_v,
                lambda_u = lambda_u,
                lambda_v = lambda_v
            )$V[, 1:rank]


            # newX should be uncencter and unscaled.
            # check new X has same colnames
            if (length(dim(newX)) != 2L) {
                moma_error("'newX' must be a matrix or data frame")
            }

            if (dim(newX)[2] != self$p) {
                moma_error(
                    paste0(
                        "`newX` is incompatible with orignal data. ",
                        "It must have ", self$p, " columns."
                    )
                )
            }

            PV <- solve(crossprod(V), t(V)) # project onto the span of V
            scaled_data <- scale(newX, self$center, self$scale)
            result <- scaled_data %*% t(PV)
            colnames(result) <- paste0("PC", seq_len(rank))

            return(list(
                scaled_data = scaled_data,
                proj_data = result,
                V = V,
                PV = PV
            ))
        }
    )
)

#' Deflation Schemes for PCA
#'
#' In \code{MoMA} three deflation schemes are provided for PCA.
#' Using terminology in the reference, they are Hotelling's deflation,
#' two-way projection deflation, and Schur complement deflation.
#'
#' See the parameter \code{deflation_scheme} argument in the function
#' \code{moma_sfpca}. Also refer to the reference below
#' for theoretical properties.
#'
#' @references Michael Weylandt. "Multi-Rank Sparse and Functional PCA: Manifold Optimization and
#' Iterative Deflation Techniques." arXiv:1907.12012v1, 2019.
#' @name PCA_deflation
NULL

#' Sparse and functional PCA
#'
#' \code{moma_sfpca} creates an \code{SFPCA} R6 object and returns it.
#' @param u_sparse,v_sparse An object of class inheriting from "\code{moma_sparsity_type}". Most conveniently
#'        specified by functions described in \code{\link{moma_sparsity_options}}. It specifies the type of sparsity-inducing
#'        penalty function used in the model. Note that for \code{moma_spca}, these two parameters must not be
#'        specified at the same time. For \code{moma_fpca} and \code{moma_twfpca}, they must not be specified.
#' @param u_smooth,v_smooth An object of class inheriting from "\code{moma_smoothness_type}". Most conveniently
#'          specified by functions described in \code{moma_smoothness}. It specifies the type of smoothness
#'           terms used in the model. Note that for \code{moma_fpca}, these two parameters must not be
#'          specified at the same time. For \code{moma_spca} and \code{moma_twspca}, they must not be specified.
#' @param deflation_scheme A string specifying the deflation scheme.
#'          It should be one of \code{"PCA_Hotelling", "PCA_Schur_Complement", "PCA_Projection"}.
#'
#' In the discussion below, let \eqn{u,v} be the normalized vectors obtained by
#' scaling the penalized singular vectors.
#'
#' When \code{deflation_scheme = "Hotelling_deflation"} is specified, the following deflation
#' scheme is used. \eqn{\boldsymbol{X}_{t} :=\boldsymbol{X}_{t-1}-d_{t} \boldsymbol{u}_{t} \boldsymbol{v}_{t}^{T}},
#' where \eqn{d_{t}=\boldsymbol{u}_{t}^{T} \boldsymbol{X}_{t-1} \boldsymbol{v}_{t}}.
#'
#' When \code{deflation_scheme = "PCA_Schur_Complement"} is specified, the following deflation
#' scheme is used: \eqn{\boldsymbol{X}_{t} :=\left(\boldsymbol{I}_{n}-
#' \boldsymbol{u}_{t} \boldsymbol{u}_{t}^{T}\right) \boldsymbol{X}_{t-1}
#' \left(\boldsymbol{I}_{p}-\boldsymbol{v}_{t} \boldsymbol{v}_{t}^{T}\right)}.
#'
#' When \code{deflation_scheme = "PCA_Projection"} is specified, the following deflation
#' scheme is used:
#' \eqn{\boldsymbol{X}_{t} :=\boldsymbol{X}_{t-1}-\frac{\boldsymbol{X}_{t-1}
#' \boldsymbol{v}_{t} \boldsymbol{u}_{t}^{T} \boldsymbol{X}_{t-1}}{\boldsymbol{u}_{t}^{T}
#' \boldsymbol{X}_{t-1} \boldsymbol{v}_{t}}}.
#' @return An R6 object which provides helper functions to access the results. See \code{\link{moma_R6}}.
#' @inheritParams moma_sfcca
#' @name moma_sfpca
#' @export
moma_sfpca <- function(X, ...,
                       center = TRUE, scale = FALSE,
                       u_sparse = moma_empty(), v_sparse = moma_lasso(),
                       u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                       pg_settings = moma_pg_settings(),
                       max_bic_iter = 5,
                       rank = 1,
                       deflation_scheme = "PCA_Hotelling") {
    chkDots(...)

    error_if_not_of_class(u_sparse, "moma_sparsity_type")
    error_if_not_of_class(v_sparse, "moma_sparsity_type")
    error_if_not_of_class(u_smooth, "moma_smoothness_type")
    error_if_not_of_class(v_smooth, "moma_smoothness_type")

    return(SFPCA$new(
        X,
        center = center, scale = scale,
        # sparsity
        u_sparsity = u_sparse$sparsity_type,
        v_sparsity = v_sparse$sparsity_type,
        lambda_u = u_sparse$lambda,
        lambda_v = v_sparse$lambda,
        # smoothness
        Omega_u = u_smooth$Omega,
        Omega_v = v_smooth$Omega,
        alpha_u = u_smooth$alpha,
        alpha_v = v_smooth$alpha,
        pg_settings = pg_settings,
        # Map strings to encoding
        select_scheme_list = list(
            select_scheme_alpha_u = SELECTION_SCHEME[[u_smooth$select_scheme]],
            select_scheme_alpha_v = SELECTION_SCHEME[[v_smooth$select_scheme]],
            select_scheme_lambda_u = SELECTION_SCHEME[[u_sparse$select_scheme]],
            select_scheme_lambda_v = SELECTION_SCHEME[[v_sparse$select_scheme]]
        ),
        max_bic_iter = max_bic_iter,
        rank = rank,
        deflation_scheme = DEFLATION_SCHEME[[deflation_scheme]]
    ))
}

#' Perform one-way sparse PCA
#'
#' \code{moma_spca} is a function for performing one-way sparse PCA.
#' @export
#' @describeIn moma_sfpca a function for performing one-way sparse PCA
moma_spca <- function(X, ...,
                      center = TRUE, scale = FALSE,
                      u_sparse = moma_empty(), v_sparse = moma_lasso(),
                      #    u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                      pg_settings = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1,
                      deflation_scheme = "PCA_Hotelling") {
    chkDots(...)
    is_u_penalized <- !missing(u_sparse)
    is_v_penalized <- !missing(v_sparse)
    if (!is_u_penalized && !is_v_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (is_u_penalized && is_v_penalized) {
        moma_error("Please use `moma_twspca` if both sides are penalized.")
    }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparse = u_sparse, v_sparse = v_sparse,
        # u_smooth = u_smooth, v_smooth = v_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank,
        deflation_scheme = deflation_scheme
    ))
    # moma_error("Not implemented: SPCA")
}


#' Perform two-way sparse PCA
#'
#' \code{moma_twspca} is a function for performing two-way sparse PCA.
#' @export
#' @describeIn moma_sfpca a function for performing two-way sparse PCA
moma_twspca <- function(X, ...,
                        center = TRUE, scale = FALSE,
                        u_sparse = moma_lasso(), v_sparse = moma_lasso(),
                        #    u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                        pg_settings = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1,
                        deflation_scheme = "PCA_Hotelling") {
    chkDots(...)
    is_u_penalized <- !missing(u_sparse)
    is_v_penalized <- !missing(v_sparse)
    if (!is_u_penalized && !is_v_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (is_u_penalized != is_v_penalized) {
        moma_warning("Please use `moma_spca` if only one side is penalized.")
    }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparse = u_sparse, v_sparse = v_sparse,
        # u_smooth = u_smooth, v_smooth = v_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank,
        deflation_scheme = deflation_scheme
    ))
}

#' Perform one-way functional PCA
#'
#' \code{moma_fpca} is a function for performing one-way functional PCA.
#' @export
#' @describeIn moma_sfpca a function for performing one-way functional PCA
moma_fpca <- function(X, ...,
                      center = TRUE, scale = FALSE,
                      #    u_sparse = moma_empty(), v_sparse = moma_empty(),
                      u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                      pg_settings = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1,
                      deflation_scheme = "PCA_Hotelling") {
    chkDots(...)
    is_u_penalized <- !missing(u_smooth)
    is_v_penalized <- !missing(v_smooth)
    if (!is_u_penalized && !is_v_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (is_u_penalized && is_v_penalized) {
        moma_error("Please use `moma_twfpca` if both sides are penalized.")
    }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparse = moma_empty(), v_sparse = moma_empty(),
        u_smooth = u_smooth, v_smooth = v_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank,
        deflation_scheme = deflation_scheme
    ))
}

#' Perform two-way functional PCA
#'
#' \code{moma_twfpca} is a function for performing two-way functional PCA.
#' @export
#' @describeIn moma_sfpca a function for performing two-way functional PCA
moma_twfpca <- function(X, ...,
                        center = TRUE, scale = FALSE,
                        #    u_sparse = moma_empty(), v_sparse = moma_empty(),
                        u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                        pg_settings = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1,
                        deflation_scheme = "PCA_Hotelling") {
    chkDots(...)
    is_u_penalized <- !missing(u_smooth)
    is_v_penalized <- !missing(v_smooth)
    if (!is_u_penalized && !is_v_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (!is_u_penalized || !is_v_penalized) {
        moma_warning("Please use `moma_fpca` if only one side is penalized.")
    }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparse = moma_empty(), v_sparse = moma_empty(),
        u_smooth = u_smooth, v_smooth = v_smooth,
        pg_settings = pg_settings,
        max_bic_iter = max_bic_iter,
        rank = rank,
        deflation_scheme = deflation_scheme
    ))
}
