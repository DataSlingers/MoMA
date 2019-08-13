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
                                      select_scheme_str = "gggg",
                                      max_bic_iter = 5,
                                      rank = 1) {
            chkDots(...)
            # Step 1: check ALL arguments
            # Step 1.1: lambdas and alphas
            if (!inherits(alpha_x, c("numeric", "integer")) ||
                !inherits(alpha_y, c("numeric", "integer")) ||
                !inherits(lambda_x, c("numeric", "integer")) ||
                !inherits(lambda_y, c("numeric", "integer"))) {
                moma_error(paste0(
                    "All penalty levels (",
                    sQuote("lambda_x"), ", ",
                    sQuote("lambda_y"), ", ",
                    sQuote("alpha_x"), ", ",
                    sQuote("alpha_y"),
                    ") must be numeric."
                ))
            }
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
            if (any(!is.finite(X))) {
                moma_error("X must not have NaN, NA, or Inf.")
            }
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
            if (!inherits(x_sparsity, "_moma_sparsity_type") || !inherits(y_sparsity, "_moma_sparsity_type")) {
                moma_error(
                    "Sparse penalty should be of class ",
                    sQuote("_moma_sparsity_type"),
                    ". Try using, for example, `x_sparsity = lasso()`."
                )
            }
            self$x_sparsity <- x_sparsity
            self$y_sparsity <- y_sparsity

            # Step 1.4: PG loop settings
            if (!inherits(pg_settings, "moma_pg_settings")) {
                moma_error(
                    "pg_settings penalty should be of class ",
                    sQuote("moma_pg_settings"),
                    ". Try using, for example, `pg_settings = moma_pg_settings(MAX_ITER=1e+4)`."
                )
            }
            self$pg_settings <- pg_settings

            # Step 1.5: smoothness
            Omega_x <- check_omega(Omega_x, alpha_x, px)
            Omega_y <- check_omega(Omega_y, alpha_y, py)
            self$Omega_x <- Omega_x
            self$Omega_y <- Omega_y

            # Step 1.6: check selection scheme string
            # "g" stands for grid search, "b" stands for BIC
            if (!inherits(select_scheme_str, "character") ||
                nchar(select_scheme_str) != 4 ||
                !all(strsplit(select_scheme_str, split = "")[[1]] %in% c("b", "g"))) {
                moma_error(
                    "Invalid select_scheme_str ", select_scheme_str,
                    ". It should be a four-char string containing only 'b' or 'g'."
                )
            }

            # turn "b"/"g" to 1/0
            # `select_scheme_list` will be passed to C++ functions
            select_scheme_list <- list(
                selection_criterion_alpha_x = 0,
                selection_criterion_alpha_y = 0,
                selection_criterion_lambda_x = 0,
                selection_criterion_lambda_y = 0
            )
            # `fixed_list` will be stored in the R6 object
            fixed_list <- list(
                # "Fixed" parameters are those
                # i) that are chosen by BIC, or
                # ii) that are not specified during initialization of the SFLDA object, or
                # iii) that are scalars as opposed to vectors during initialization of the SFLDA object.
                is_alpha_x_fixed = FALSE,
                is_alpha_y_fixed = FALSE,
                is_lambda_x_fixed = FALSE,
                is_lambda_y_fixed = FALSE
            )

            parameter_length_list <- sapply(FUN = length, list(
                self$alpha_x,
                self$alpha_y,
                self$lambda_x,
                self$lambda_y
            ))
            for (i in 1:4) {
                select_scheme_list[[i]] <- ifelse(substr(select_scheme_str, i, i) == "g", 0, 1)
                fixed_list[[i]] <- substr(select_scheme_str, i, i) == "b" || parameter_length_list[i] == 1
            }
            self$select_scheme_list <- select_scheme_list
            self$fixed_list <- fixed_list

            # Step 1.7: check rank
            # TODO: check that `rank` < min(rank(X), rank(Y))
            # w.r.t to certain numric precision
            if (!inherits(rank, "numeric") ||
                !is.wholenumber(rank) ||
                rank <= 0 ||
                rank > min(px, py, n)) { # LDA_SPECIAL_PART
                moma_error("`rank` should be a positive integer smaller than the rank of the data matrix.")
            }
            self$rank <- rank

            # Step 2: pack all arguments in a list
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
                    selection_criterion_alpha_u = select_scheme_list$selection_criterion_alpha_x,
                    selection_criterion_alpha_v = select_scheme_list$selection_criterion_alpha_y,
                    selection_criterion_lambda_u = select_scheme_list$selection_criterion_lambda_x,
                    selection_criterion_lambda_v = select_scheme_list$selection_criterion_lambda_y
                ),
                list(
                    max_bic_iter = max_bic_iter
                ),
                list(
                    deflation_scheme = DEFLATION_SCHEME["LDA"] # LDA_SPECIAL_PART
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
            parameter_length_list <- sapply(FUN = length, list(
                alpha_x,
                alpha_y,
                lambda_x,
                lambda_y
            ))
            if (any(parameter_length_list > 1)) {
                moma_error("Non-length-one input in SFLDA::get_mat_by_index.")
            }

            # indices should be integers
            if (!all(
                is.wholenumber(alpha_x),
                is.wholenumber(alpha_y),
                is.wholenumber(lambda_x),
                is.wholenumber(lambda_y)
            )) {
                moma_error("Non-integer input in SFLDA::get_mat_by_index.")
            }

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
                is_fixed <- self$fixed_list
                if (any(is_fixed == TRUE & is_missing == FALSE)) {
                    moma_error(
                        paste0(
                            "Invalid index in SFLDA::get_mat_by_index. Do not specify indexes of parameters ",
                            "i) that are chosen by BIC, or ",
                            "ii) that are not specified during initialization of the SFLDA object, or ",
                            "iii) that are scalars during initialization of the SFLDA object."
                        )
                    )
                }
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
                is_fixed <- self$fixed_list
                if (any(is_fixed == TRUE & is_missing == FALSE)) {
                    moma_error(
                        paste0(
                            "Invalid index in SFLDA::left_project. Do not specify indexes of parameters ",
                            "i) that are chosen by BIC, or ",
                            "ii) that are not specified during initialization of the SFLDA object, or ",
                            "iii) that are scalars during initialization of the SFLDA object."
                        )
                    )
                }
            }

            if (rank > self$rank) {
                moma_error("Invalid `rank` in SFLDA::left_project.")
            }

            parameter_length_list <- sapply(FUN = length, list(
                alpha_x,
                alpha_y,
                lambda_x,
                lambda_y
            ))
            if (any(parameter_length_list > 1) ||
                !all(
                    is.wholenumber(alpha_x),
                    is.wholenumber(alpha_y),
                    is.wholenumber(lambda_x),
                    is.wholenumber(lambda_y)
                )) {
                moma_error("Non-integer input in SFLDA::left_project.")
            }

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
        },

        plot = function() {
            shinyApp(
                ui = fluidPage(
                    tags$style(
                        type = "text/css",
                        ".recalculating { opacity: 1.0; }"
                    ),
                    titlePanel("Sparse and functional LDA"),
                    sidebarLayout(
                        sidebarPanel(
                            width = 2,
                            sliderInput(
                                "alpha_x_i",
                                "alpha_x",
                                min = 1,
                                max = ifelse(self$fixed_list$is_alpha_x_fixed,
                                    1, length(self$alpha_x)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "alpha_y_i",
                                "alpha_y",
                                min = 1,
                                max = ifelse(self$fixed_list$is_alpha_y_fixed,
                                    1, length(self$alpha_y)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "lambda_x_i",
                                "lambda_x",
                                min = 1,
                                max = ifelse(self$fixed_list$is_lambda_x_fixed,
                                    1, length(self$lambda_x)
                                ),
                                value = 1,
                                step = 1
                            ),
                            sliderInput(
                                "lambda_y_i",
                                "lambda_y",
                                min = 1,
                                max = ifelse(self$fixed_list$is_lambda_y_fixed,
                                    1, length(self$lambda_y)
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
                                        column(width = 5, plotOutput("X_PC_loadings_plot")),
                                        column(width = 5, plotOutput("Y_group_scores_plot"))
                                    ),
                                    fluidRow(
                                        column(4, verbatimTextOutput("alpha_x")),
                                        column(4, verbatimTextOutput("lambda_x")),
                                        column(4, verbatimTextOutput("alpha_y")),
                                        column(4, verbatimTextOutput("lambda_y"))
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
                            alpha_x = input$alpha_x_i,
                            alpha_y = input$alpha_y_i,
                            lambda_x = input$lambda_x_i,
                            lambda_y = input$lambda_y_i
                        )
                    })

                    get_left_projected <- reactive({
                        private$private_X_project(
                            newX = self$X,
                            alpha_x = input$alpha_x_i,
                            alpha_y = input$alpha_y_i,
                            lambda_x = input$lambda_x_i,
                            lambda_y = input$lambda_y_i,
                            rank = 2
                        )
                    })

                    # plot the singuler vector
                    output$Y_group_scores_plot <- renderPlot({
                        k <- as.integer(input$rank)
                        rank_k_result <- get_rank_k_result()
                        plot(rank_k_result$Y_group_scores[, k],
                            ylab = "Y_group_scores",
                            type = "l",
                            xaxt = "n",
                            main = "Scores of groups",
                            xlab = "Groups"
                        )
                        axis(1,
                            at = 1:self$py,
                            labels = self$y_coln
                        )
                    })

                    # plot the singuler vector
                    output$X_PC_loadings_plot <- renderPlot({
                        k <- as.integer(input$rank)
                        rank_k_result <- get_rank_k_result()
                        plot(rank_k_result$X_PC_loadings[, k],
                            ylab = "X_PC_loadings",
                            type = "l",
                            xaxt = "n",
                            main = "Plot of X loadings",
                            xlab = "Features of X"
                        )
                        axis(1,
                            at = 1:self$px,
                            labels = self$x_coln
                        )
                    })

                    output$lambda_y <- renderPrint({
                        k <- as.integer(input$rank)

                        alpha_x_value <- get_rank_k_result()$chosen_alpha_x[k]
                        cat(paste0("alpha_x = ", alpha_x_value), "\n")
                        alpha_y_value <- get_rank_k_result()$chosen_alpha_y[k]
                        cat(paste0("alpha_y = ", alpha_y_value), "\n")
                        lambda_x_value <- get_rank_k_result()$chosen_lambda_x[k]
                        cat(paste0("lambda_x = ", lambda_x_value), "\n")
                        lambda_y_value <- get_rank_k_result()$chosen_lambda_y[k]
                        cat(paste0("lambda_y = ", lambda_y_value), "\n")
                    })

                    output$X_rows_projected <- renderPlot({
                        X_projected <- get_left_projected()$proj_data
                        # print(X_projected)
                        plot(X_projected,
                            xlab = "PC1",
                            ylab = "PC2",
                            main = "Projected rows of the data matrix",
                            cex = 0
                        )
                        text(X_projected, labels = self$Y_factor, cex = 0.7)
                    })

                    output$Y_rows_projected <- renderPlot({

                    })
                }
            )
        }
    )
)


#' Perform two-way sparse and functional LDA
#'
#' \code{moma_sflda} creates an \code{SFLDA} R6 object and returns.
#' @param X,Y_factor data matrix.
#' @param ... force users to specify arguments by names
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' Defaults to \code{TRUE}.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance.
#' Defaults to \code{FALSE}.
#' @param x_sparse,y_sparse an object of class inheriting from "\code{moma_sparsity_type}". Most conveniently
#'        specified by functions described in \code{\link{moma_sparsity}}. It specifies the type of sparsity-inducing
#'        penalty function used in the model. Note that for \code{moma_slda}, these two parameter must not be
#'        specified at the same time. For \code{moma_flda} and \code{moma_twflda}, they must not be specified.
#' @param x_smooth,y_smooth an object of class inheriting from "\code{moma_smoothness_type}". Most conveniently
#'          specified by functions described in \code{moma_smoothness}. It specifies the type of smoothness
#'           terms used in the model. Note that for \code{moma_flda}, these two parameter must not be
#'          specified at the same time. For \code{moma_slda} and \code{moma_twslda}, they must not be specified.
#' @param pg_setting an object of class inheriting from "\code{moma_sparsity}". Most conviently
#'          specified by functions described in \code{\link{moma_pg_settings}}. It specifies the type of algorithm
#'          used to solve the problem, acceptable level of precision, and the maximum number of iterations allowed.
#' @param max_bic_iter a positive integer. Defaults to 5. The maximum number of iterations allowed
#' in nested greedy BIC selection scheme.
#' @param rank a positive integer. Defaults to 1. The maximal rank, i.e., maximal number of principal components to be used.
#' @export
moma_sflda <- function(X, ..., Y_factor,
                       center = TRUE, scale = FALSE,
                       x_sparse = moma_empty(), y_sparse = moma_empty(),
                       x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                       pg_settings = moma_pg_settings(),
                       max_bic_iter = 5,
                       rank = 1) {
    chkDots(...)
    if (!inherits(x_sparse, "moma_sparsity_type") ||
        !inherits(y_sparse, "moma_sparsity_type")) {
        moma_error(
            "Invalid argument: ",
            sQuote("x_sparse / y_sparse"),
            ". They should be of class `moma_sparsity_type`."
        )
    }

    if (!inherits(x_smooth, "moma_smoothness_type") ||
        !inherits(y_smooth, "moma_smoothness_type")) {
        moma_error(
            "Invalid argument: ",
            sQuote("x_smooth / y_smooth"),
            ". They should be of class `moma_smoothness_type`."
        )
    }


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
        select_scheme_str = paste0( # the order is important
            x_smooth$select_scheme,
            y_smooth$select_scheme,
            x_sparse$select_scheme,
            y_sparse$select_scheme
        ),
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}

#' Perform one-way sparse LDA
#'
#' \code{moma_slda} is a wrapper around R6 object \code{SFLDA}
#' @export
#' @describeIn moma_sflda a function for one-way sparse LDA
moma_slda <- function(X, ..., Y_factor,
                      center = TRUE, scale = FALSE,
                      x_sparse = moma_empty(), y_sparse = moma_empty(),
                      #    x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                      pg_setting = moma_pg_settings(),
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
#' \code{moma_twslda} is a wrapper around R6 object \code{SFLDA}
#' @export
#' @describeIn moma_sflda a function for two-way sparse LDA
moma_twslda <- function(X, ..., Y_factor,
                        center = TRUE, scale = FALSE,
                        x_sparse = moma_empty(), y_sparse = moma_empty(),
                        #    x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                        pg_setting = moma_pg_settings(),
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
#' \code{moma_flda} is a wrapper around R6 object \code{SFLDA}
#' @export
#' @describeIn moma_sflda a function for one-way functional LDA
moma_flda <- function(X, ..., Y_factor,
                      center = TRUE, scale = FALSE,
                      #    x_sparse = moma_empty(), y_sparse = moma_empty(),
                      x_smooth = moma_smoothness(), y_smooth = moma_smoothness(),
                      pg_setting = moma_pg_settings(),
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
#' \code{moma_twflda} is a wrapper around R6 object \code{SFLDA}
#' @export
#' @describeIn moma_sflda a function for two-way functional LDA
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
