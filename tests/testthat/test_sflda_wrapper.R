context("Test SFLDA R interface")

test_that("SFLDA object: Error when X and Y are incompatible", {
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- factor(sample(num_group, size = n - 1, replace = TRUE)) # error
    expect_error(
        SFLDA$new(X = X, Y_factor = Y),
        "`X` and `Y_factor` must have the same number of samples."
    )
})

test_that("SFLDA object: Error when Y is not a factor", {
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10

    # `Y` is not a factor
    Y <- sample(num_group, size = n, replace = TRUE)
    expect_error(
        SFLDA$new(X = X, Y_factor = Y),
        "`Y` must be a factor"
    )

    expect_no_error(
        SFLDA$new(X = X, Y_factor = as.factor(Y))
    )

    # wrong argument name: it should be `Y_factor`
    expect_warning(
        expect_error(
            SFLDA$new(X = X, Y = as.factor(Y)),
            "argument \"Y_factor\" is missing, with no default"
        ),
        paste0(sQuote("Y"), " will be disregarded")
    )
})

test_that("SFLDA object: `scale_X` default to FALSE", {
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- as.factor(sample(num_group, size = n, replace = TRUE))

    expect_equal(
        SFLDA$new(X = X, Y_factor = Y)$scale_X,
        FALSE
    )
})

test_that("SFLDA object: Consistent with direct SVD result", {
    set.seed(12)
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- as.factor(sample(num_group, size = n, replace = TRUE))

    # run LDA with SFLDA$new
    a <- SFLDA$new(X = X, Y_factor = Y, center = FALSE)$grid_result[[1]]
    au <- a$u$vector
    av <- a$v$vector
    expect_equal(norm(au, "F"), 1)
    expect_equal(norm(av, "F"), 1)

    # run LDA with SVD
    scatter_mat <- t(X) %*% model.matrix(~ Y - 1) / sqrt(n)
    scatter_mat_svd <- svd(scatter_mat)

    # compare
    expect_equal(
        au,
        matrix(scatter_mat_svd$u[, 1])
    )
    expect_equal(
        av,
        matrix(scatter_mat_svd$v[, 1])
    )
    expect_equal(
        a$d,
        scatter_mat_svd$d[1]
    )
})

test_that("SFLDA object: Correct deflation scheme", {
    set.seed(12)
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- as.factor(sample(num_group, size = n, replace = TRUE))

    a <- SFLDA$new(X = X, Y_factor = Y, center = FALSE, rank = 3)$grid_result

    scatter_mat <- t(X) %*% model.matrix(~ Y - 1) / sqrt(n)

    rank1 <- get_5Dlist_elem(a, 1, 1, 1, 1, 1)[[1]]
    rank2 <- get_5Dlist_elem(a, 1, 1, 1, 1, 2)[[1]]
    rank3 <- get_5Dlist_elem(a, 1, 1, 1, 1, 3)[[1]]


    # This function deflates `lda_list$X`
    # by the outer product of `lda_list$u$vector`
    find_deflated_mat <- function(lda_list) {
        x_scores <- lda_list$u$vector
        # canonical covariates
        x_cv <- lda_list$X %*% x_scores
        x_cv <- x_cv / norm(x_cv, "F")

        return(lda_list$X - x_cv %*% t(x_cv) %*% lda_list$X)
    }

    # This function checks that
    # `lda_list$u$vector` and `lda_list$v$vector`
    # are constructed by finding pSVD of the scatter
    # matrix `lda_list$X %*% Y`.
    check_lda_uv <- function(lda_list) {
        au <- lda_list$u$vector
        av <- lda_list$v$vector
        X <- lda_list$X

        expect_equal(norm(au, "F"), 1)
        expect_equal(norm(av, "F"), 1)

        # Y is defined outside this function.
        scatter_mat <- t(X) %*% model.matrix(~ Y - 1) / sqrt(n)
        scatter_mat_svd <- svd(scatter_mat)

        expect_equal(
            au,
            matrix(scatter_mat_svd$u[, 1])
        )
        expect_equal(
            av,
            matrix(scatter_mat_svd$v[, 1])
        )
        expect_equal(
            lda_list$d,
            scatter_mat_svd$d[1]
        )
    }

    # check that `rank1`, `rank2` and
    # `rank3` are self-consistent.
    check_lda_uv(rank1)
    check_lda_uv(rank2)
    check_lda_uv(rank3)

    # check that `rank2$X` and `rank1$X`,
    # `rank2$X` and `rank3$X` are consistent
    expect_equal(
        rank1$X,
        X
    )
    expect_equal(
        rank2$X,
        find_deflated_mat(rank1)
    )
    expect_equal(
        rank3$X,
        find_deflated_mat(rank2)
    )
})

test_that("SFLDA object: Compare with the built-in MASS::lda", {

})

test_that("SFLDA object: Column / row names of a named matrix is stored correctly", {
    px <- 4
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))


    Y <- factor(c("a", "a", "b", "b", "c", "c"))
    a <- SFLDA$new(X = X, Y_factor = Y)

    expect_equal(
        a$x_coln,
        colnames(X)
    )
    expect_equal(
        a$x_rown,
        rownames(X)
    )
    expect_equal(
        a$y_coln,
        levels(Y)
    )
})

test_that("SFLDA object: get_mat_by_id", {
    px <- 4
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))


    Y <- factor(c("a", "a", "b", "b", "c", "c"))
    a <- SFLDA$new(
        X = X, Y_factor = Y, rank = 2,
        alpha_x = c(1, 2, 3)
    )

    expect_no_error(
        a$get_mat_by_index()
    )
    expect_no_error(
        a$get_mat_by_index(alpha_x = 2)
    )
    expect_no_error(
        a$get_mat_by_index(alpha_x = 3)
    )

    # check `get_mat_by_index` defaults to
    # returning `get_mat_by_index(1,1,1,1)`
    expect_equal(
        a$get_mat_by_index(),
        a$get_mat_by_index(alpha_x = 1)
    )

    # check the `chosen_alpha_x` is correct
    for (i in 1:3)
    {
        expect_equal(
            a$get_mat_by_index(alpha_x = i)$chosen_alpha_x,
            c(i, i)
        )
    }

    # check the list returned by `get_mat_by_index`
    # contains expected contents
    expect_true(
        all(
            c(
                "X_PC_loadings", "Y_group_scores", "d", "chosen_lambda_x",
                "chosen_lambda_y", "chosen_alpha_x", "chosen_alpha_y"
            ) %in%
                names(a$get_mat_by_index())
        )
    )
})

test_that("SFLDA object: print", {
    px <- 4
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))


    Y <- factor(c("a", "a", "b", "b", "c", "c"))
    a <- SFLDA$new(
        X = X, Y_factor = Y, rank = 2,
        alpha_x = c(1, 2, 3)
    )

    print_message <- capture.output(print(a))

    expected_message <- c(
        "An <SFLDA> object containing solutions to the following settings",
        "Rank:  2 ",
        "Penalty and selection:",
        "alpha_x: grid search, range: 1 2 3 ",
        "alpha_y: grid search, range: 0 ",
        "lambda_x: grid search, range: 0 ",
        "lambda_y: grid search, range: 0 "
    )
    expect_equal(expected_message, print_message)
})

test_that("SFLDA object: select_scheme", {
    px <- 4
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))


    Y <- factor(c("a", "a", "b", "b", "c", "c"))
    a <- SFLDA$new(
        X = X, Y_factor = Y, rank = 2,
        alpha_x = c(1, 2, 3),
        alpha_y = seq(0, 1, 0.4),
        lambda_x = seq(0, 1, 0.41),
        lambda_y = seq(0, 1, 0.42)
    )

    # Omega's default to second diff mat.
    expect_equal(a$Omega_x, second_diff_mat(4))
    expect_equal(a$Omega_y, second_diff_mat(3))

    # check the default selection scheme is grid search
    expect_true(all(
        a$select_scheme_list == c(0, 0, 0, 0)
    ))


    a <- SFLDA$new(
        X = X, Y_factor = Y, rank = 2,
        alpha_x = c(1, 2, 3),
        alpha_y = seq(0, 1, 0.4),
        lambda_x = seq(0, 1, 0.41),
        lambda_y = seq(0, 1, 0.42),
        select_scheme_list = list(
            select_scheme_alpha_x = SELECTION_SCHEME[["bic"]],
            select_scheme_alpha_y = SELECTION_SCHEME[["grid"]],
            select_scheme_lambda_x = SELECTION_SCHEME[["bic"]],
            select_scheme_lambda_y = SELECTION_SCHEME[["grid"]]
        )
    )

    # check the selection scheme has been correctly
    # specified
    expect_true(all(
        a$select_scheme_list == c(1, 0, 1, 0)
    ))


    a <- SFLDA$new(
        X = X, Y_factor = Y, rank = 2,
        alpha_x = c(1, 2, 3),
        alpha_y = seq(0, 1, 0.4),
        lambda_x = seq(0, 1, 0.41),
        lambda_y = seq(0, 1, 0.42),
        select_scheme_list = list(
            select_scheme_alpha_x = SELECTION_SCHEME[["bic"]],
            select_scheme_alpha_y = SELECTION_SCHEME[["bic"]],
            select_scheme_lambda_x = SELECTION_SCHEME[["bic"]],
            select_scheme_lambda_y = SELECTION_SCHEME[["bic"]]
        )
    )

    # check the selection scheme has been correctly
    # specified
    expect_true(all(
        a$select_scheme_list == c(1, 1, 1, 1)
    ))
})


test_that("SFLDA special-case functions", {
    px <- 4
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))


    Y <- factor(c("a", "a", "b", "b", "c", "c"))
    expect_no_error(
        a <- moma_sflda(X,
            Y_factor = Y,
            x_sparse = moma_lasso(lambda = seq(0, 1, 0.4)),
            y_smooth = moma_smoothness(select_scheme = "b")
        )
    )

    expect_true(all(
        a$select_scheme_list == c(0, 1, 0, 0)
    ))
})
