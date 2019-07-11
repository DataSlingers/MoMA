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
#'   \item{\code{get_mat}}{Arguments: \code{alpha_u}, \code{alpha_v}, \code{lambda_u}, \code{lambda_v}
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
SFPCA <- R6::R6Class("SFPCA", list(
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
    n = NULL,
    p = NULL,
    X = NULL,
    initialize = function(X, ...,
                              center = TRUE, scale = FALSE,
                              u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                              Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                              pg_setting = moma_pg_settings(),
                              selection_scheme_str = "gggg",
                              max_bic_iter = 5,
                              rank = 1) {

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
        n <- dim(X)[1]
        p <- dim(X)[2]
        X <- scale(X, center = center, scale = scale)

        cen <- attr(X, "scaled:center")
        sc <- attr(X, "scaled:scale")
        if (any(sc == 0)) {
            moma_error("cannot rescale a constant/zero column to unit variance")
        }

        self$center <- if (is.null(cen)) FALSE else cen
        self$scale <- if (is.null(sc)) FALSE else sc
        self$n <- n
        self$p <- p
        self$X <- X

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

        # Step 1.5: smoothness
        Omega_u <- check_omega(Omega_u, alpha_u, n)
        Omega_v <- check_omega(Omega_v, alpha_v, p)
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
        for (i in 1:4) {
            selection_scheme_list[[i]] <- ifelse(strsplit(selection_scheme_str, split = "")[[1]][i] == "g", 0, 1)
        }
        self$selection_scheme_list <- selection_scheme_list

        # Step 1.7: check rank
        if (!inherits(rank, "numeric") || rank <= 0
        || rank > min(p, n)) {
            moma_error("rank should be a legit positive integer.")
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
                Omega_u = check_omega(Omega_u, alpha_u, n),
                Omega_v = check_omega(Omega_v, alpha_v, p),
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

    get_mat = function(alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1) {
        # Sanity check: if a parameter has been chosen by BIC, then
        # the index for that parameter should not be specified.
        parameter <- c(alpha_u, alpha_v, lambda_u, lambda_v)
        if (any(parameter != 1 & self$selection_scheme_list != 0)) {
            moma_error("Invalid index in SFPCA::get_mat. Do not specify indexes of parameters chosen by BIC.")
        }


        n <- self$n
        p <- self$p
        rank <- self$rank

        U <- matrix(0, nrow = n, ncol = rank)
        V <- matrix(0, nrow = p, ncol = rank)
        for (i in (1:self$rank)) {
            U[, i] <-
                get_5Dlist_elem(self$grid_result,
                    alpha_u_i = alpha_u,
                    lambda_u_i = lambda_u,
                    alpha_v_i = alpha_v,
                    lambda_v_i = lambda_v, rank_i = i
                )[[1]]$u$vector

            V[, i] <-
                get_5Dlist_elem(self$grid_result,
                    alpha_u_i = alpha_u,
                    lambda_u_i = lambda_u,
                    alpha_v_i = alpha_v,
                    lambda_v_i = lambda_v, rank_i = i
                )[[1]]$v$vector
        }

        coln <- if (is.null(colnames(self$X))) paste0("Xcol_", seq_len(p)) else colnames(self$X)
        rown <- if (is.null(rownames(self$X))) paste0("Xrow_", seq_len(n)) else rownames(self$X)
        dimnames(V) <-
            list(coln, paste0("PC", seq_len(rank)))
        dimnames(U) <-
            list(rown, paste0("PC", seq_len(rank)))
        return(list(U = U, V = V))
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
        cat(paste0("alpha_u: ", selection_list_str[2]), "\n")
        print(self$alpha_v)
        cat(paste0("lambda_u: ", selection_list_str[3]), "\n")
        print(self$lambda_u)
        cat(paste0("lambda_v: ", selection_list_str[4]), "\n")
        print(self$lambda_v)
    },

    left_project = function(newX, ...,
                                alpha_u = 1, alpha_v = 1, lambda_u = 1, lambda_v = 1) {

        # check indexes
        parameter <- c(alpha_u, alpha_v, lambda_u, lambda_v)
        if (any(parameter != 1 & self$selection_scheme_list != 0)) {
            moma_error("Invalid index in SFPCA::get_mat. Do not specify indexes of parameters chosen by BIC.")
        }
        V <- self$get_mat(
            alpha_u = alpha_u,
            alpha_v = alpha_v,
            lambda_u = lambda_u,
            lambda_v = lambda_v
        )$V


        # newX should be uncencter and unscaled.
        # check new X has same colnames
        if (length(dim(newX)) != 2L) {
            moma_error("'newX' must be a matrix or data frame")
        }

        if (dim(newX)[2] != self$p) {
            moma_error("`newX` is incompatible with orignal data.")
        }

        PV <- solve(crossprod(V), t(V))
        scaled_data <- scale(newX, self$center, self$scale)
        result <- scaled_data %*% t(PV)
        colnames(result) <- paste0("PC", seq_len(self$rank))
        return(list(
            scaled_data = scaled_data,
            proj_data = result,
            V = V,
            PV = PV
        ))
    }
))

#' Perform two-way sparse and functional PCA
#'
#' \code{moma_sfpca} creates an \code{SFPCA} R6 object and returns.
#' @param X data matrix.
#' @param ... force users to specify arguments by names
#' @param center a logical value indicating whether the variables should be shifted to be zero centered.
#' Defaults to \code{TRUE}.
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance.
#' Defaults to \code{FALSE}.
#' @param u_sparsity,v_sparsity an object of class inheriting from "\code{moma_sparsity}". Most conveniently
#' specified by functions described in \code{\link{moma_sparsity}}. It specifies the type of sparsity-inducing
#' penalty function used in the model. Note that for \code{moma_spca}, these two parameter must not be
#' specified at the same time. For \code{moma_fpca} and \code{moma_twfpca}, they must not be specified.
#' @param lambda_u,lambda_v a numeric vector or a number that specifies the penalty level of sparsity.
#' Note that for \code{moma_spca}, these two parameter must not be
#' specified at the same time. For \code{moma_fpca} and \code{moma_twfpca}, they must not be specified.
#' @param Omega_u,Omega_v a positive definite matrix that encourages smoothness.  Note that for \code{moma_fpca}, these two parameter must not be
#' specified at the same time. For \code{moma_spca} and \code{moma_twspca}, they must not be specified.
#' @param alpha_u,alpha_v v a numeric vector or a number that specifies the penalty level of smoothness.
#' Note that for \code{moma_fpca}, these two parameter must not be
#' specified at the same time. For \code{moma_spca} and \code{moma_twspca}, they must not be specified.
#' @param pg_setting an object of class inheriting from "\code{moma_sparsity}". Most conviently
#' specified by functions described in \code{\link{moma_pg_settings}}. It specifies the type of algorithm
#' used to solve the problem, acceptable level of precision, and the maximum number of iterations allowed.
#' @param selection_scheme_str a one-letter, two-letter or four-letter string that specifies selection schemes for tuning
#' parameters, containing only "b" and "g". "b" stands for greedy nested BIC selection, and "g"
#' stands for exhaustive grid search. For \code{moma_sfpca}, it is a four-letter string, and selection schemes for the tuning parameters,
#' alpha_u, alpha_v, lambda_u and lambda_v are specified by the four letters of the string, respectively. For
#' \code{moma_spca} and \code{moma_fpca}, it is a one-leter string that specifies the selection scheme
#' for the paramter of interest. For \code{moma_twspca}, it is a two-letter string, and
#' selection schemes for lambda_u and lambda_v are specified by the two letters respectively. For
#' \code{moma_twfpca}, it is a two-letter string, and selection schemes for alpha_u and alpha_v
#' are specified by the two letters respectively.
#' @param max_bic_iter a positive integer. Defaults to 5. The maximum number of iterations allowed
#' in nested greedy BIC selection scheme.
#' @param rank a positive integer. Defaults to 1. The maximal rank, i.e., maximal number of principal components to be used.
#' @export
moma_sfpca <- function(X, ...,
                       center = TRUE, scale = FALSE,
                       u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                       Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0, # so is alpha_u/_v
                       pg_setting = moma_pg_settings(),
                       selection_scheme_str = "gggg",
                       max_bic_iter = 5,
                       rank = 1) {
    chkDots(...)
    return(SFPCA$new(
        X,
        center = center, scale = scale,
        # sparsity
        u_sparsity = u_sparsity, v_sparsity = v_sparsity,
        lambda_u = lambda_u, lambda_v = lambda_v,
        # smoothness
        Omega_u = Omega_u, Omega_v = Omega_v,
        alpha_u = alpha_u, alpha_v = alpha_v,
        pg_setting = pg_setting,
        selection_scheme_str = selection_scheme_str,
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
                      u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                      pg_setting = moma_pg_settings(),
                      selection_scheme_str = "g",
                      max_bic_iter = 5,
                      rank = 1) {
    chkDots(...)
    u_penalized <- !(missing(u_sparsity) && missing(lambda_u))
    v_penalized <- !(missing(v_sparsity) && missing(lambda_v))
    if (!u_penalized && !v_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (u_penalized && v_penalized) {
        moma_error("Please use `moma_twspca` if both sides are penalized.")
    }

    if (nchar(selection_scheme_str) != 1 || !selection_scheme_str %in% c("g", "b")) {
        moma_error("`selection_scheme_str` should be either 'g' or 'b'")
    }

    # selection_scheme_str is in order of alpha_u/v, lambda_u/v
    full_selection_scheme_str <-
        if (u_penalized) {
            paste0("gg", selection_scheme_str, "g")
        }
        else if (v_penalized) {
            paste0("ggg", selection_scheme_str)
        }
        else {
            "gggg"
        }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparsity = u_sparsity, v_sparsity = v_sparsity,
        lambda_u = lambda_u, lambda_v = lambda_v,
        #  Omega_u = Omega_u, Omega_v = Omega_v,
        #  alpha_u = alpha_u, alpha_v = alpha_v,
        pg_setting = pg_setting,
        selection_scheme_str = full_selection_scheme_str,
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
                        u_sparsity = empty(), v_sparsity = empty(), lambda_u = 0, lambda_v = 0, # lambda_u/_v is a vector or scalar
                        pg_setting = moma_pg_settings(),
                        selection_scheme_str = "gg",
                        max_bic_iter = 5,
                        rank = 1) {
    chkDots(...)
    u_penalized <- !(missing(u_sparsity) && missing(lambda_u))
    v_penalized <- !(missing(v_sparsity) && missing(lambda_v))
    if (!u_penalized && !v_penalized) {
        moma_warning("No sparsity is imposed!")
    }

    if (!u_penalized || !v_penalized) {
        moma_warning("Please use `moma_spca` if only one side is penalized.")
    }

    if (nchar(selection_scheme_str) != 2) {
        moma_error("`selection_scheme_str` should be of length two.")
    }
    if (!all(strsplit(selection_scheme_str, split = "")[[1]] %in% c("b", "g"))) {
        moma_error("`selection_scheme_str` should consist of 'g' or 'b'")
    }

    # selection_scheme_str is in order of alpha_u/v, lambda_u/v
    full_selection_scheme_str <- paste0("gg", selection_scheme_str)

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        u_sparsity = u_sparsity, v_sparsity = v_sparsity,
        lambda_u = lambda_u, lambda_v = lambda_v,
        # Omega_u = Omega_u, Omega_v = Omega_v,
        # alpha_u = alpha_u, alpha_v = alpha_v,
        pg_setting = pg_setting,
        selection_scheme_str = full_selection_scheme_str,
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
                      Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0,
                      pg_setting = moma_pg_settings(),
                      selection_scheme_str = "g",
                      max_bic_iter = 5,
                      rank = 1) {
    chkDots(...)
    u_penalized <- !(missing(Omega_u) && missing(alpha_u))
    v_penalized <- !(missing(Omega_v) && missing(alpha_v))
    if (!u_penalized && !v_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (u_penalized && v_penalized) {
        moma_error("Please use `moma_twfpca` if both sides are penalized.")
    }

    if (nchar(selection_scheme_str) != 1 || !selection_scheme_str %in% c("g", "b")) {
        moma_error("`selection_scheme_str` should be either 'g' or 'b'")
    }

    # selection_scheme_str is in order of alpha_u/v, lambda_u/v
    full_selection_scheme_str <-
        if (u_penalized) {
            paste0(selection_scheme_str, "ggg")
        }
        else if (v_penalized) {
            paste0("g", selection_scheme_str, "gg")
        }
        else {
            "gggg"
        }

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        #   u_sparsity = u_sparsity, v_sparsity = v_sparsity,
        #   lambda_u = lambda_u, lambda_v = lambda_v,
        Omega_u = Omega_u, Omega_v = Omega_v,
        alpha_u = alpha_u, alpha_v = alpha_v,
        pg_setting = pg_setting,
        selection_scheme_str = full_selection_scheme_str,
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
                        Omega_u = NULL, Omega_v = NULL, alpha_u = 0, alpha_v = 0,
                        pg_setting = moma_pg_settings(),
                        selection_scheme_str = "gg",
                        max_bic_iter = 5,
                        rank = 1) {
    chkDots(...)
    u_penalized <- !(missing(Omega_u) && missing(alpha_u))
    v_penalized <- !(missing(Omega_v) && missing(alpha_v))
    if (!u_penalized && !v_penalized) {
        moma_warning("No smoothness is imposed!")
    }

    if (!u_penalized || !v_penalized) {
        moma_warning("Please use `moma_fpca` if only one side is penalized.")
    }

    if (nchar(selection_scheme_str) != 2) {
        moma_error("`selection_scheme_str` should be of length two.")
    }
    if (!all(strsplit(selection_scheme_str, split = "")[[1]] %in% c("b", "g"))) {
        moma_error("`selection_scheme_str` should consist of 'g' or 'b'")
    }

    # selection_scheme_str is in order of alpha_u/v, lambda_u/v
    full_selection_scheme_str <- paste0(selection_scheme_str, "gg")

    return(moma_sfpca(
        X,
        center = center, scale = scale,
        #   u_sparsity = u_sparsity, v_sparsity = v_sparsity,
        #   lambda_u = lambda_u, lambda_v = lambda_v,
        Omega_u = Omega_u, Omega_v = Omega_v,
        alpha_u = alpha_u, alpha_v = alpha_v,
        pg_setting = pg_setting,
        selection_scheme_str = full_selection_scheme_str,
        max_bic_iter = max_bic_iter,
        rank = rank
    ))
}
