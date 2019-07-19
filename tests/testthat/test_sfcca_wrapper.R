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
        alpha_u = 1, alpha_v = 1,
        Omega_u = second_diff_mat(px),
        Omega_v = second_diff_mat(py),
        lambda_u = 1, lambda_v = 1,
        u_sparsity = lasso(), v_sparsity = scad()
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
        x_cv <- cca_list$X %*% x_scores / norm(cca_list$X %*% x_scores, "F")
        y_cv <- cca_list$Y %*% y_scores / norm(cca_list$Y %*% y_scores, "F")
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
    alpha_u <- seq(0, 1, 0.21)
    alpha_v <- seq(0, 1, 0.32)
    lambda_u <- seq(0, 1, 0.43)
    lambda_v <- seq(0, 1, 0.54)

    SFCCA$undebug("initialize")
    a <- SFCCA$new(
        X = X, Y = Y,
        center = FALSE, scale = FALSE,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(px),
        Omega_v = second_diff_mat(py),
        lambda_u = lambda_u, lambda_v = lambda_v,
        u_sparsity = lasso(), v_sparsity = scad()
    )

    for (av in 1:length(alpha_v)) {
        for (au in 1:length(alpha_u)) {
            for (lv in 1:length(lambda_v)) {
                for (lu in 1:length(lambda_u)) {
                    single_call <-
                        SFCCA$new(
                            X = X, Y = Y,
                            center = FALSE, scale = FALSE,
                            alpha_u = alpha_u[au],
                            alpha_v = alpha_v[av],
                            Omega_u = second_diff_mat(px),
                            Omega_v = second_diff_mat(py),
                            lambda_u = lambda_u[lu],
                            lambda_v = lambda_v[lv],
                            u_sparsity = lasso(),
                            v_sparsity = scad()
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
