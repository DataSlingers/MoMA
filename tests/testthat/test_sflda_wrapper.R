context("Test SFLDA R interface")

test_that("Error when X and Y are incompatible", {
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

test_that("Error when Y is not a factor", {
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- sample(num_group, size = n, replace = TRUE)
    expect_error(
        SFLDA$new(X = X, Y_factor = Y),
        "`Y` must be a factor"
    )

    expect_no_error(
        SFLDA$new(X = X, Y_factor = as.factor(Y))
    )

    # it should be `Y_factor`
    expect_warning(
        expect_error(
            SFLDA$new(X = X, Y = as.factor(Y)),
            "argument \"Y_factor\" is missing, with no default"
        ),
        paste0(
            "extra argument ",
            sQuote("Y"),
            " will be disregarded"
        )
    )
})

test_that("`scale_X` default to FALSE", {
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

test_that("Consistent with direct SVD result", {
    set.seed(12)
    px <- 4
    num_group <- 5 # number of groups
    n <- 6
    X <- matrix(runif(n * px), n, px) * 10
    Y <- as.factor(sample(num_group, size = n, replace = TRUE))

    a <- SFLDA$new(X = X, Y_factor = Y, center = FALSE)$grid_result[[1]]
    au <- a$u$vector
    av <- a$v$vector
    expect_equal(norm(au, "F"), 1)
    expect_equal(norm(av, "F"), 1)

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
        a$d,
        scatter_mat_svd$d[1]
    )
})

test_that("Correct deflation scheme", {
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

    find_deflated_mat <- function(lda_list) {
        x_scores <- lda_list$u$vector
        # canonical covariates
        x_cv <- lda_list$X %*% x_scores / norm(lda_list$X %*% x_scores, "F")
        return(lda_list$X - x_cv %*% t(x_cv) %*% lda_list$X)
    }

    # This function checks that
    # u and v are constructed by
    # finding pSVD of the scatter
    # matrix, which is constructed
    # by
    check_lda_uv <- function(lda_list) {
        au <- lda_list$u$vector
        av <- lda_list$v$vector
        X <- lda_list$X

        expect_equal(norm(au, "F"), 1)
        expect_equal(norm(av, "F"), 1)

        # Y is defined outside this function.
        # It is unchanged
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

    check_lda_uv(rank1)
    check_lda_uv(rank2)
    check_lda_uv(rank3)

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

test_that("Compare with the built-in MASS::lda", {

})
