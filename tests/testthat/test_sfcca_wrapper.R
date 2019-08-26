context("Test SFCCA R interface")

test_that("Solve a penalized CCA with three canonical covariates", {
    set.seed(1)
    px <- 4
    py <- 5
    n <- 100
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10

    a <- SFCCA$new(
        X = X, Y = Y, rank = 3,
        center = FALSE, scale = FALSE,
        alpha_x = 1, alpha_y = 1,
        Omega_x = second_diff_mat(px),
        Omega_y = second_diff_mat(py),
        lambda_x = 1, lambda_y = 1,
        x_sparsity = lasso(), y_sparsity = scad()
    )

    # contains X and Y
    grid_names <- names(a$grid_result[[1]])
    if (
        !("X" %in% grid_names) ||
            !("Y" %in% grid_names)
    ) {
        moma_error("`SFCCA` returns an unexpected list.")
    }

    rank1 <- get_5Dlist_elem(a$grid_result, 1, 1, 1, 1, 1)[[1]]
    rank2 <- get_5Dlist_elem(a$grid_result, 1, 1, 1, 1, 2)[[1]]
    rank3 <- get_5Dlist_elem(a$grid_result, 1, 1, 1, 1, 3)[[1]]

    find_deflated_mat <- function(cca_list) {
        x_scores <- cca_list$u$vector
        y_scores <- cca_list$v$vector
        # canonical covariates
        x_cv <- cca_list$X %*% x_scores
        x_cv <- x_cv / norm(x_cv, "F")
        y_cv <- cca_list$Y %*% y_scores
        y_cv <- y_cv / norm(y_cv, "F")

        list(
            X = cca_list$X - x_cv %*% t(x_cv) %*% cca_list$X,
            Y = cca_list$Y - y_cv %*% t(y_cv) %*% cca_list$Y
        )
    }

    expect_equal(
        rank1$X,
        X
    )
    expect_equal(
        rank1$Y,
        Y
    )
    expect_equal(
        rank2$X,
        find_deflated_mat(rank1)$X
    )
    expect_equal(
        rank2$Y,
        find_deflated_mat(rank1)$Y
    )
    expect_equal(
        rank3$X,
        find_deflated_mat(rank2)$X
    )
    expect_equal(
        rank3$Y,
        find_deflated_mat(rank2)$Y
    )
})

test_that("Evaluate on a grid", {
    set.seed(1)
    px <- 4
    py <- 5
    n <- 100
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10
    alpha_x <- seq(0, 1, 0.21)
    alpha_y <- seq(0, 1, 0.32)
    lambda_x <- seq(0, 1, 0.43)
    lambda_y <- seq(0, 1, 0.54)

    SFCCA$undebug("initialize")
    a <- SFCCA$new(
        X = X, Y = Y,
        center = FALSE, scale = FALSE,
        alpha_x = alpha_x, alpha_y = alpha_y,
        Omega_x = second_diff_mat(px),
        Omega_y = second_diff_mat(py),
        lambda_x = lambda_x, lambda_y = lambda_y,
        x_sparsity = lasso(), y_sparsity = scad()
    )

    for (av in 1:length(alpha_y)) {
        for (au in 1:length(alpha_x)) {
            for (lv in 1:length(lambda_y)) {
                for (lu in 1:length(lambda_x)) {
                    single_call <-
                        SFCCA$new(
                            X = X, Y = Y,
                            center = FALSE, scale = FALSE,
                            alpha_x = alpha_x[au],
                            alpha_y = alpha_y[av],
                            Omega_x = second_diff_mat(px),
                            Omega_y = second_diff_mat(py),
                            lambda_x = lambda_x[lu],
                            lambda_y = lambda_y[lv],
                            x_sparsity = lasso(),
                            y_sparsity = scad()
                        )$grid_result[[1]]

                    target <- get_5Dlist_elem(
                        a$grid_result,
                        au, lu, av, lv, 1
                    )[[1]]
                    expect_equal(single_call, target)
                }
            }
        }
    }
})

test_that("CCA as LDA, orthogonal design", {

    # Reference:
    # Penalized Discriminant Analysis,
    # Trevor Hastie, Andreas Buja, and Robert Tibshirani,
    # Section 3.2 Penalized canonical correlation analysis

    set.seed(12)
    px <- 4
    py <- 5
    num_group <- 5 # number of groups
    n <- 30

    # test for 20 data sets
    for (i in 1:20) {
        # Note for non-orthogonal design
        # such equivalance does not establish
        # X <- matrix(runif(n * px), n, px)
        # Y <- matrix(runif(n * py), n, py)

        # generate orthogonal designs
        X <- eigen(crossprod(matrix(runif(n * n), n, n)))$vector[, 1:px]
        Y <- eigen(crossprod(matrix(runif(n * n), n, n)))$vector[, 1:py]

        # form within-class and between-class variance
        sig11 <- 1 / n * t(Y) %*% Y
        sig22 <- 1 / n * t(X) %*% X
        sig12 <- 1 / n * t(Y) %*% X
        sig21 <- t(sig12)
        M <- solve(sig11, sig12) # whose rows are means
        sig_bet <- t(M) %*% sig11 %*% M
        sig_w <- sig22 - sig_bet
        L <- chol(sig_w) # L %*% t(L) = sig_w
        L_inv <- solve(L)

        # first discriminant covariate
        first_dc <- L_inv %*% eigen(t(L_inv) %*% sig_bet %*% L_inv)$vector[, 1]
        first_dc <- first_dc / norm(first_dc, "F")

        # my first right canonical covariate
        my_first <- SFCCA$new(X = X, Y = Y, center = FALSE, scale = FALSE)$grid_result[[1]]$u$vector

        expect_equal(
            abs(first_dc),
            abs(my_first)
        )
    }
})

test_that("SFCCA object: Error when X and Y are incompatible", {
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif((n + 1) * py), n + 1, py) * 10

    expect_error(
        SFCCA$new(X = X, Y = Y),
        "`X` and `Y` must have the same number of samples."
    )
})


test_that("SFCCA object: Consistent with direct SVD result", {
    set.seed(12)
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10

    # run CCA with SFCCA$new
    a <- SFCCA$new(X = X, Y = Y, center = FALSE)$grid_result[[1]]
    au <- a$u$vector
    av <- a$v$vector
    expect_equal(norm(au, "F"), 1)
    expect_equal(norm(av, "F"), 1)

    # run CCA with SVD
    scatter_mat <- t(X) %*% Y
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


test_that("SFCCA object: Correct deflation scheme", {
    set.seed(12)
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10

    a <- SFCCA$new(X = X, Y = Y, center = FALSE, rank = 3)$grid_result

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

        y_scores <- lda_list$v$vector
        # canonical covariates
        y_cv <- lda_list$Y %*% y_scores
        y_cv <- y_cv / norm(y_cv, "F")

        newX <- lda_list$X - x_cv %*% t(x_cv) %*% lda_list$X
        newY <- lda_list$Y - y_cv %*% t(y_cv) %*% lda_list$Y
        return(list(
            newX, newY
        ))
    }

    # This function checks that
    # `lda_list$u$vector` and `lda_list$v$vector`
    # are constructed by finding pSVD of the scatter
    # matrix `lda_list$X %*% Y`.
    check_lda_uv <- function(lda_list) {
        au <- lda_list$u$vector
        av <- lda_list$v$vector
        X <- lda_list$X
        Y <- lda_list$Y

        expect_equal(norm(au, "F"), 1)
        expect_equal(norm(av, "F"), 1)

        scatter_mat <- t(X) %*% Y
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
        rank1$Y,
        Y
    )
    expect_equal(
        rank2$X,
        find_deflated_mat(rank1)[[1]]
    )
    expect_equal(
        rank2$Y,
        find_deflated_mat(rank1)[[2]]
    )
    expect_equal(
        rank3$X,
        find_deflated_mat(rank2)[[1]]
    )
    expect_equal(
        rank3$Y,
        find_deflated_mat(rank2)[[2]]
    )
})


test_that("SFCCA object: Column / row names of a named matrix is stored correctly", {
    set.seed(123)
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10
    colnames(X) <- paste0("XFeature_", seq(px))
    rownames(X) <- paste0("sample_", seq(n))

    a <- SFCCA$new(X = X, Y = Y)

    # column names of X
    expect_equal(
        a$x_coln,
        colnames(X)
    )
    expect_equal(
        a$x_rown,
        rownames(X)
    )

    # column names of Y are set to defualt values
    expect_equal(
        a$y_coln,
        paste0("Ycol_", seq(1, py))
    )
})


test_that("SFLDA object: get_mat_by_id", {
    set.seed(123)
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10

    a <- SFCCA$new(
        X = X, Y = Y, rank = 2,
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
                "X_PC_loadings", "Y_PC_loadings", "d", "chosen_lambda_x",
                "chosen_lambda_y", "chosen_alpha_x", "chosen_alpha_y"
            ) %in%
                names(a$get_mat_by_index())
        )
    )
})


test_that("SFCCA object: print", {
    set.seed(123)
    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10

    a <- SFCCA$new(
        X = X, Y = Y, rank = 2,
        alpha_x = c(1, 2, 3)
    )

    print_message <- capture.output(print(a))

    expected_message <- c(
        "An <SFCCA> object containing solutions to the following settings",
        "Rank:  2 ",
        "Penalty and selection:",
        "alpha_x: grid search, range: 1 2 3 ",
        "alpha_y: grid search, range: 0 ",
        "lambda_x: grid search, range: 0 ",
        "lambda_y: grid search, range: 0 "
    )
    expect_equal(expected_message, print_message)
})

test_that("SFCCA special-case functions: get_mat_by_index", {
    set.seed(123)

    px <- 4
    py <- 5
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- matrix(runif(n * py), n, py) * 10
    a <- moma_sfcca(
        X = X, Y = Y,
        x_sparse = moma_lasso(lambda = seq(0, 2, 0.11), select_scheme = "b")
    )

    expect_true(all(a$fixed_list == c(1, 1, 1, 1)))

    expect_no_error(a$get_mat_by_index())
    expect_error(a$get_mat_by_index(lambda_x = 1))

    expect_no_error(a$X_project(X))
    expect_no_error(a$Y_project(Y))

    expect_error(a$Y_project(X))
    expect_error(a$X_project(Y))
})
