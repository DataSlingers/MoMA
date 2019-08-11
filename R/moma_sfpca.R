#' \code{SFPCA} \code{R6} object
#'
#' An \code{R6} object for performing SFPCA and parameter selection
#'
#' During initialziation of an \code{SFPCA} object, \code{R}
#' calls the \code{C++}-side function, \code{cpp_multirank_BIC_grid_search}, and
#' wraps the results returned. The \code{SFPCA} object also records penalty levels
#' and selection schemes of tuning parameters. Several helper
#' methods are provivded to facilitate access to results.
#' Initialization is delegated to \code{\link{moma_sfpca}}.
#' @seealso \code{\link{moma_sfpca}},
#' \code{\link{moma_spca}},
#' \code{\link{moma_fpca}},
#' \code{\link{moma_twspca}},
#' \code{\link{moma_twfpca}}
#'
#' @section Members:
#'
#' \describe{
#'   \item{\code{center,scale}}{The attributes "\code{scaled:center}" and "\code{scaled:scale}" of function \code{scale}.
#' The numeric centering and scalings used (if any) of the data matrix.}
#'   \item{\code{grid_result}}{a 5-D list containing the results evaluated on the paramter grid.}
#'   \item{\code{selection_scheme_list}}{a list with elements \code{selection_criterion_alpha_u},
#'            \code{selection_criterion_alpha_v},
#'            \code{selection_criterion_lambda_u},
#'            \code{selection_criterion_lambda_v}. Each of them is either 0 or 1. 0 stands for grid search
#' and 1 stands for BIC search. TODO: descibe the mixed selection procedure.}
#' }
#' @section Methods:
#'
#' \describe{
#'   \item{\code{get_mat_by_index}}{Arguments: \code{alpha_u}, \code{alpha_v}, \code{lambda_u}, \code{lambda_v}
#' are the indices of the parameters in the paramter grid, which is specified during initialization.
#'
#' Obtain the right and left sigular penalized vectors, which are packed into matrices \code{U} and \code{V}.}
#'   \item{\code{print}}{Display tuning parameters and selection schemes.}
#'   \item{\code{left_project}}{Arguments: \code{newX}, a matrix (un-centered and un-scaled) of
#' the same number of columns as the original data matrix.
#'
#' Project the new data into the space spaned by the
#' penalized left singular vectors, after scaling and centering as needed.}
#' }
#' @return an \code{SFPCA} R6 object.
#' @export
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
        selection_scheme_list = NULL,
        pg_setting = NULL,
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
                                      pg_setting = moma_pg_settings(),
                                      selection_scheme_str = "gggg",
                                      max_bic_iter = 5,
                                      rank = 1) {
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

            # Step 1.3: sparsity
            error_if_not_of_class(u_sparsity, "moma_sparsity")
            error_if_not_of_class(v_sparsity, "moma_sparsity")
            self$u_sparsity <- u_sparsity
            self$v_sparsity <- v_sparsity

            # Step 1.4: PG loop settings
            error_if_not_of_class(pg_setting, "moma_pg_settings")
            self$pg_setting <- pg_setting

            # Step 1.5: smoothness
            Omega_u <- check_omega(Omega_u, alpha_u, n)
            Omega_v <- check_omega(Omega_v, alpha_v, p)
            self$Omega_u <- Omega_u
            self$Omega_v <- Omega_v

            # Step 1.6: check selection scheme string
            # "g" stands for grid search, "b" stands for BIC
            error_if_not_fourchar_bg_string(selection_scheme_str)

            # turn "b"/"g" to 1/0
            # `selection_scheme_list` will be passed to C++ functions
            selection_scheme_list <- list(
                selection_criterion_alpha_u = 0,
                selection_criterion_alpha_v = 0,
                selection_criterion_lambda_u = 0,
                selection_criterion_lambda_v = 0
            )
            # `fixed_list` will be stored in the R6 object
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
            parameter_length_list <- vapply(FUN = length, list(
                self$alpha_u,
                self$alpha_v,
                self$lambda_u,
                self$lambda_v
            ), integer(1))

            for (i in 1:4) {
                para_select_str_i <- substr(selection_scheme_str, i, i)
                selection_scheme_list[[i]] <- ifelse(para_select_str_i == "g", 0, 1)
                fixed_list[[i]] <- para_select_str_i == "b" || parameter_length_list[[i]] == 1
            }
            self$selection_scheme_list <- selection_scheme_list
            self$fixed_list <- fixed_list

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
                pg_setting,
                selection_scheme_list,
                list(
                    max_bic_iter = max_bic_iter
                )
            )
            # make sure we explicitly specify ALL arguments
            if (length(algo_settings_list) != length(formals(cpp_multirank_BIC_grid_search))) {
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

            error_if_not_finite_numeric_scalar(alpha_u)
            error_if_not_finite_numeric_scalar(alpha_v)
            error_if_not_finite_numeric_scalar(lambda_u)
            error_if_not_finite_numeric_scalar(lambda_v)

            error_if_not_wholenumber(alpha_u)
            error_if_not_wholenumber(alpha_v)
            error_if_not_wholenumber(lambda_u)
            error_if_not_wholenumber(lambda_v)

            # A "fixed" parameter should not be specified
            # at all (this is a bit stringent, can be improved later).
            # "Fixed" parameters are those
            # i) that are chosen by BIC, or
            # ii) that are not specified during initialization of the SFCPA object, or
            # iii) that are scalars as opposed to vectors during initialization of the SFCPA object.

            # When `get_mat_by_index` is called internally
            # we skip the input checking
            if (private$check_input_index) {
                is_missing <- list(missing(alpha_u), missing(alpha_v), missing(lambda_u), missing(lambda_v))
                is_fixed <- self$fixed_list
                if (any(is_fixed == TRUE & is_missing == FALSE)) {
                    moma_error(
                        paste0(
                            "Invalid index in SFPCA::get_mat_by_index. Do not specify indexes of parameters ",
                            "i) that are chosen by BIC, or ",
                            "ii) that are not specified during initialization of the SFCPA object, or ",
                            "iii) that are scalars during initialization of the SFCPA object."
                        )
                    )
                }
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
            if (any(self$selection_scheme_list != 0)) {
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
            is_fixed <- self$fixed_list == TRUE
            param_str_list <- c("alpha_u", "alpha_v", "lambda_u", "lambda_v")
            if (any(is_missing == FALSE & is_fixed == TRUE)) { # elementwise and
                output_para <- is_missing == FALSE & is_fixed == TRUE
                moma_error(
                    paste0(
                        "Invalid index in SFPCA::interpolate: ",
                        paste(param_str_list[output_para], collapse = ", "),
                        ". Do not specify indexes of parameters ",
                        "i) that are chosen by BIC, or ",
                        "ii) that are not specified during initialization of the SFCPA object, or ",
                        "iii) that are scalars during initialization of the SFCPA object."
                    )
                )
            }
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
                    pg_setting = self$pg_setting,
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
            selection_list_str <- lapply(self$selection_scheme_list, function(x) {
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
                is_missing <- list(missing(alpha_u), missing(alpha_v), missing(lambda_u), missing(lambda_v))
                is_fixed <- self$fixed_list
                if (any(is_fixed == TRUE & is_missing == FALSE)) {
                    moma_error(
                        paste0(
                            "Invalid index in SFPCA::get_mat_by_index. Do not specify indexes of parameters ",
                            "i) that are chosen by BIC, or ",
                            "ii) that are not specified during initialization of the SFCPA object, or ",
                            "iii) that are scalars during initialization of the SFCPA object."
                        )
                    )
                }
            }

            if (rank > self$rank) {
                moma_error("Invalid `rank` in SFPCA::left_project.")
            }

            error_if_not_finite_numeric_scalar(alpha_u)
            error_if_not_finite_numeric_scalar(alpha_v)
            error_if_not_finite_numeric_scalar(lambda_u)
            error_if_not_finite_numeric_scalar(lambda_v)

            error_if_not_wholenumber(alpha_u)
            error_if_not_wholenumber(alpha_v)
            error_if_not_wholenumber(lambda_u)
            error_if_not_wholenumber(lambda_v)

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
        },

        plot = function() {
            shinyApp(
                ui = fluidPage(
                    tags$style(
                        type = "text/css",
                        ".recalculating { opacity: 1.0; }"
                    ),
                    titlePanel("Sparse and functional PCA"),
                    sidebarLayout(
                        sidebarPanel(
                            width = 2,
                            sliderInput(
                                "alpha_u_i",
                                "alpha_u",
                                min = 1,
                                max = ifelse(self$fixed_list$is_alpha_u_fixed,
                                    1, length(self$alpha_u)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "alpha_v_i",
                                "alpha_v",
                                min = 1,
                                max = ifelse(self$fixed_list$is_alpha_v_fixed,
                                    1, length(self$alpha_v)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "lambda_u_i",
                                "lambda_u",
                                min = 1,
                                max = ifelse(self$fixed_list$is_lambda_u_fixed,
                                    1, length(self$lambda_u)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "lambda_v_i",
                                "lambda_v",
                                min = 1,
                                max = ifelse(self$fixed_list$is_lambda_v_fixed,
                                    1, length(self$lambda_v)
                                ),
                                value = 1,
                                step = 1
                            )
                        ),
                        mainPanel(
                            tabsetPanel(
                                tabPanel(
                                    "Loadings of PCs",
                                    fluidRow(
                                        column(
                                            width = 2,
                                            radioButtons("rank", "Rank", seq(1, self$rank))
                                        ),
                                        column(width = 5, plotOutput("u_loadings_plot")),
                                        column(width = 5, plotOutput("v_loadings_plot"))
                                    ),
                                    fluidRow(
                                        column(4, verbatimTextOutput("alpha_u")),
                                        column(4, verbatimTextOutput("lambda_u")),
                                        column(4, verbatimTextOutput("alpha_v")),
                                        column(4, verbatimTextOutput("lambda_v"))
                                    )
                                ),
                                tabPanel(
                                    "Projected data",
                                    fluidRow(
                                        column(width = 6, plotOutput("X_rows_projected")),
                                        column(width = 6, plotOutput("X_cols_projected"))
                                    )
                                )
                            )
                        )
                    )
                ),
                server = function(input, output) {
                    get_rank_k_result <- reactive({
                        private$private_get_mat_by_index(
                            alpha_u = input$alpha_u_i,
                            alpha_v = input$alpha_v_i,
                            lambda_u = input$lambda_u_i,
                            lambda_v = input$lambda_v_i
                        )
                    })

                    get_left_projected <- reactive({
                        private$private_left_project(
                            newX = self$X,
                            alpha_u = input$alpha_u_i,
                            alpha_v = input$alpha_v_i,
                            lambda_u = input$lambda_u_i,
                            lambda_v = input$lambda_v_i,
                            rank = 2
                        )
                    })

                    # plot the singuler vector
                    output$v_loadings_plot <- renderPlot({
                        k <- as.integer(input$rank)
                        rank_k_result <- get_rank_k_result()
                        plot(rank_k_result$V[, k],
                            ylab = "v",
                            type = "l",
                            xaxt = "n"
                        )
                        axis(1,
                            at = 1:self$p,
                            labels = self$X_coln
                        )
                    })

                    # plot the singuler vector
                    output$u_loadings_plot <- renderPlot({
                        k <- as.integer(input$rank)
                        rank_k_result <- get_rank_k_result()
                        plot(rank_k_result$U[, k],
                            ylab = "u",
                            type = "l"
                        )
                    })

                    output$lambda_v <- renderPrint({
                        k <- as.integer(input$rank)

                        alpha_u_value <- get_rank_k_result()$chosen_alpha_u[k]
                        cat(paste0("alpha_u = ", alpha_u_value), "\n")
                        alpha_v_value <- get_rank_k_result()$chosen_alpha_v[k]
                        cat(paste0("alpha_v = ", alpha_v_value), "\n")
                        lambda_u_value <- get_rank_k_result()$chosen_lambda_u[k]
                        cat(paste0("lambda_u = ", lambda_u_value), "\n")
                        lambda_v_value <- get_rank_k_result()$chosen_lambda_v[k]
                        cat(paste0("lambda_v = ", lambda_v_value), "\n")
                    })

                    output$X_rows_projected <- renderPlot({
                        X_projected <- get_left_projected()$proj_data
                        # print(X_projected)
                        plot(X_projected,
                            xlab = "PC1",
                            ylab = "PC2",
                            main = "Projected rows of the data matrix"
                        )
                    })

                    output$X_cols_projected <- renderPlot({

                    })
                }
            )
        }
    )
)

#' Perform two-way sparse and functional PCA
#'
#' \code{moma_sfpca} creates an \code{SFPCA} R6 object and returns.
#' @param X data matrix.
#' @param ... force users to specify arguments by names
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' Defaults to \code{TRUE}.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance.
#' Defaults to \code{FALSE}.
#' @param u_sparse,v_sparse an object of class inheriting from "\code{moma_sparsity_type}". Most conveniently
#'        specified by functions described in \code{\link{moma_sparsity}}. It specifies the type of sparsity-inducing
#'        penalty function used in the model. Note that for \code{moma_spca}, these two parameter must not be
#'        specified at the same time. For \code{moma_fpca} and \code{moma_twfpca}, they must not be specified.
#' @param u_smooth,v_smooth an object of class inheriting from "\code{moma_smoothness_type}". Most conveniently
#'          specified by functions described in \code{moma_smoothness}. It specifies the type of smoothness
#           terms used in the model. Note that for \code{moma_fpca}, these two parameter must not be
#'          specified at the same time. For \code{moma_spca} and \code{moma_twspca}, they must not be specified.
#' @param pg_setting an object of class inheriting from "\code{moma_sparsity}". Most conviently
#'          specified by functions described in \code{\link{moma_pg_settings}}. It specifies the type of algorithm
#'          used to solve the problem, acceptable level of precision, and the maximum number of iterations allowed.
#' @param max_bic_iter a positive integer. Defaults to 5. The maximum number of iterations allowed
#' in nested greedy BIC selection scheme.
#' @param rank a positive integer. Defaults to 1. The maximal rank, i.e., maximal number of principal components to be used.
#' @export
moma_sfpca <- function(X, ...,
                       center = TRUE, scale = FALSE,
                       u_sparse = moma_empty(), v_sparse = moma_lasso(),
                       u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                       pg_setting = moma_pg_settings(),
                       max_bic_iter = 5,
                       rank = 1) {
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
        pg_setting = pg_setting,
        selection_scheme_str = paste0( # the order is important
            u_smooth$select_scheme,
            v_smooth$select_scheme,
            u_sparse$select_scheme,
            v_sparse$select_scheme
        ),
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform one-way sparse PCA
#'
#' \code{moma_spca} is a wrapper around R6 object \code{SFPCA}
#' @export
#' @describeIn moma_sfpca a function for one-way sparse PCA
moma_spca <- function(X, ...,
                      center = TRUE, scale = FALSE,
                      u_sparse = moma_empty(), v_sparse = moma_lasso(),
                      #    u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                      pg_setting = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1) {
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
        pg_setting = pg_setting,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
    # moma_error("Not implemented: SPCA")
}


#' Perform two-way sparse PCA
#'
#' \code{moma_twspca} is a wrapper around R6 object \code{SFPCA}
#' @export
#' @describeIn moma_sfpca a function for two-way sparse PCA
moma_twspca <- function(X, ...,
                        center = TRUE, scale = FALSE,
                        u_sparse = moma_lasso(), v_sparse = moma_lasso(),
                        #    u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                        pg_setting = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1) {
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
        pg_setting = pg_setting,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform one-way functional PCA
#'
#' \code{moma_fpca} is a wrapper around R6 object \code{SFPCA}
#' @export
#' @describeIn moma_sfpca a function for one-way functional PCA
moma_fpca <- function(X, ...,
                      center = TRUE, scale = FALSE,
                      #    u_sparse = moma_empty(), v_sparse = moma_empty(),
                      u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                      pg_setting = moma_pg_settings(),
                      max_bic_iter = 5,
                      rank = 1) {
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
        pg_setting = pg_setting,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform two-way functional PCA
#'
#' \code{moma_twfpca} is a wrapper around R6 object \code{SFPCA}
#' @export
#' @describeIn moma_sfpca a function for two-way functional PCA
moma_twfpca <- function(X, ...,
                        center = TRUE, scale = FALSE,
                        #    u_sparse = moma_empty(), v_sparse = moma_empty(),
                        u_smooth = moma_smoothness(), v_smooth = moma_smoothness(),
                        pg_setting = moma_pg_settings(),
                        max_bic_iter = 5,
                        rank = 1) {
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
        pg_setting = pg_setting,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}
