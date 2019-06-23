context("Test L1 linear trend filtering")
test_that("Return original data when lambda = 0", {
    set.seed(33)
    rep <- 23
    p <- 29
    for (i in 1:rep) {
        x <- runif(p)
        expect_equal(test_prox_l1gf(x, 0), matrix(t(x)))
    }
})

test_that("Return linear regression fit when lambda = infty", {
    set.seed(33)
    rep <- 23
    p <- 29
    infty <- 1e+5
    for (i in 1:rep) {
        x <- runif(p)
        t <- 1:p
        lr_fit <- predict(lm(x ~ t, data.frame(t = t)))
        expect_equal(test_prox_l1gf(x, infty), matrix(lr_fit))
    }
})

sec_diff_mat <- function(m) {
    D <- matrix(rep(0, m * (m - 2)), m - 2)
    for (i in 1:m - 2) {
        D[i, i] <- 1
        D[i, (i + 1)] <- -2
        D[i, (i + 2)] <- 1
    }
    return(D)
}

test_that("Compare with `CVX` with random lambda's", {
    p <- 10

    # answers obtained from matlab package `CVX`
    x <- c(.22, .11, .11, .06, .4, .45, .37, .76, .63, .77)
    # WARNING: probably it it not accurate itself
    ans <- matrix(c(
        0.220000, 0.176667, 0.151667, 0.126667, 0.109894, 0.095532, 0.081170, 0.066809, 0.052447, 0.038364, 0.038364,
        0.110000, 0.146667, 0.146667, 0.146667, 0.142553, 0.137234, 0.131915, 0.126596, 0.121277, 0.116061, 0.116061,
        0.110000, 0.116667, 0.141667, 0.166667, 0.175213, 0.178936, 0.182660, 0.186383, 0.190106, 0.193758, 0.193758,
        0.060000, 0.160000, 0.190143, 0.196786, 0.207872, 0.220638, 0.233404, 0.246170, 0.258936, 0.271455, 0.271455,
        0.400000, 0.286404, 0.291000, 0.295000, 0.302979, 0.312249, 0.321520, 0.330790, 0.340061, 0.349152, 0.349152,
        0.450000, 0.399298, 0.391857, 0.393214, 0.398085, 0.403860, 0.409635, 0.415410, 0.421185, 0.426848, 0.426848,
        0.370000, 0.512193, 0.492714, 0.491429, 0.493191, 0.495471, 0.497751, 0.500030, 0.502310, 0.504545, 0.504545,
        0.760000, 0.625088, 0.593571, 0.589643, 0.588298, 0.587082, 0.585866, 0.584650, 0.583435, 0.582242, 0.582242,
        0.630000, 0.694035, 0.691429, 0.687857, 0.683404, 0.678693, 0.673982, 0.669271, 0.664559, 0.659939, 0.659939,
        0.770000, 0.762982, 0.789286, 0.786071, 0.778511, 0.770304, 0.762097, 0.753891, 0.745684, 0.737636, 0.737636
    ),
    ncol = 10
    )
    lambda_set <- seq(0, 0.5, 0.05)

    cnt <- 1
    for (lambda in lambda_set) {
        my_l1tf <- round(test_prox_l1gf(x, lambda), 2)

        # At least the first two digits are the same.
        expect_lt(sum((my_l1tf - as.matrix(ans[cnt, ]))^2), 0.6e-3)

        cnt <- cnt + 1
    }
})


test_that("Commutability with affine adjustment", {
    # Ref:l1 Trend Filtering by Seung-Jean Kim
    set.seed(33)

    rep <- 23
    p <- 29
    infty <- 1e+5

    alpha <- runif(1)
    beta <- runif(1)

    for (i in 1:rep) {
        x <- runif(p)

        # remove a line from x
        t <- 1:p
        x_ <- x - alpha * t - beta
        for (lambda in seq(0, 3, 0.2)) {
            expect_equal(
                test_prox_l1gf(x_, lambda),
                test_prox_l1gf(x, lambda) - alpha * t - beta
            )
        }
    }
})

test_that("Eqivalent to fused lasso when using first-diff-mat", {
    set.seed(2)
    rep <- 5
    p <- 17
    for (i in 1:rep) {
        x <- runif(p)
        for (lambda in seq(0, 5, 0.3)) {
            # WARNING: primal-dual method only attains
            # low precision

            # primal-dual methods
            pd <- test_prox_l1gf(x, lambda, 0)

            # path algorithm
            pa <- test_prox_fusedlassopath(x, lambda)

            # Most of the results are the same
            # for at least 5 digits after decimal points,
            # but some give outrageous errors
            expect_lte(sum((pd - pa)^2), 1e-9)
        }
    }
})

test_that("Tests for difference matrix", {
    mat1 <- matrix(c(
        1, -1, 0, 0,
        0, 1, -1, 0,
        0, 0, 1, -1
    ), byrow = TRUE, nrow = 3)

    mat2 <- matrix(c(
        -1, 2, -1, 0, 0,
        0, -1, 2, -1, 0,
        0, 0, -1, 2, -1
    ), byrow = TRUE, nrow = 3)

    mat3 <- matrix(c(
        1, -5, 10, -10, 5, -1, 0, 0, 0, 0,
        0, 1, -5, 10, -10, 5, -1, 0, 0, 0,
        0, 0, 1, -5, 10, -10, 5, -1, 0, 0,
        0, 0, 0, 1, -5, 10, -10, 5, -1, 0,
        0, 0, 0, 0, 1, -5, 10, -10, 5, -1
    ), byrow = TRUE, nrow = 5)

    mat4 <- matrix(c(
        -1, 4, -6, 4, -1, 0, 0, 0, 0, 0,
        0, -1, 4, -6, 4, -1, 0, 0, 0, 0,
        0, 0, -1, 4, -6, 4, -1, 0, 0, 0,
        0, 0, 0, -1, 4, -6, 4, -1, 0, 0,
        0, 0, 0, 0, -1, 4, -6, 4, -1, 0,
        0, 0, 0, 0, 0, -1, 4, -6, 4, -1
    ), byrow = TRUE, nrow = 6)

    mat5 <- matrix(c(
        1, -3, 3, -1, 0, 0, 0,
        0, 1, -3, 3, -1, 0, 0,
        0, 0, 1, -3, 3, -1, 0,
        0, 0, 0, 1, -3, 3, -1
    ), byrow = TRUE, nrow = 4)
    expect_equal(norm(mat1 - l1tf_diff_mat(4, 0)), 0)
    expect_equal(norm(mat2 - l1tf_diff_mat(5, 1)), 0)
    expect_equal(norm(mat5 - l1tf_diff_mat(7, 2)), 0)
    expect_equal(norm(mat4 - l1tf_diff_mat(10, 3)), 0)
    expect_equal(norm(mat3 - l1tf_diff_mat(10, 4)), 0)
})

test_that("Equivalent to polynomial regression", {
    # NOTE: theory to be confirmed

    set.seed(22)
    large_lam <- 100000

    t <- seq(0, 10, 0.1)

    # an "N"-shape curve
    y <- (t - 2)^2 + 0.2 * (3 - t)^3 + rnorm(length(t), sd = 2)

    for (deg in seq(2)) {
        # WARNING: tests fail for deg >= 5
        # polynomial regression fit
        md <- lm(y ~ poly(t, deg, raw = TRUE))
        pr <- rep(0, length(t))
        for (i in 1:(deg + 1)) {
            # poly regression
            pr <- pr + coef(md)[i] * t^(i - 1)
        }

        # trend filtering fit
        tf <- test_prox_l1gf(y, 1000000, deg)

        expect_equal(matrix(pr), tf)
    }
})
