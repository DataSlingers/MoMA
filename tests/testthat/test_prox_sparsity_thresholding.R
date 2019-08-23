context("Thresolding Tests")
test_that("Non-convexity parameter checks", {
    lambda <- 1
    x <- seq(-4, 4, 0.01)

    # Check valid non-convexity parameter for SCAD, MCP

    # Throw error if below threshold
    expect_error(test_prox_scad(x, lambda, 1.99))
    expect_error(test_prox_mcp(x, lambda, 0.99))
    expect_error(test_prox_scad(x, lambda, 2))
    expect_error(test_prox_mcp(x, lambda, 1))
    # Otherwise succeed
    expect_no_error(test_prox_scad(x, lambda, 2.1))
    expect_no_error(test_prox_mcp(x, lambda, 1.1))
})

test_that("When lambda = 0, prox operators are no-ops", {
    lambda <- 0
    x <- matrix(seq(-4, 4, 0.01), ncol = 1)

    for (prox_func in c(
        test_prox_lasso,
        test_prox_scad,
        test_prox_mcp,
        test_prox_scadvec,
        test_prox_mcpvec,
        test_prox_fusedlassopath
    )) {
        expect_equal(
            x,
            prox_func(x, lambda)
        )
    }
})

test_that("When lambda = 0, non-negative prox operators zero-out negative values", {
    lambda <- 0
    x <- matrix(seq(-4, 4, 0.01), ncol = 1)

    for (prox_func in c(
        test_prox_nnlasso,
        test_prox_nnscad,
        test_prox_nnmcp
    )) {
        expect_equal(
            x * (x >= 0),
            prox_func(x, lambda)
        )
    }
})

test_that("Prox operators return correct results for lambda = 3", {
    lambda <- 1
    gamma <- 3
    x <- seq(-4, 4, 0.05)

    # Worked out by hand, gamma = 3 for non-convex operators
    lasso.goal <- matrix(c(-3.00, -2.95, -2.90, -2.85, -2.80, -2.75, -2.70, -2.65, -2.60, -2.55, -2.50, -2.45, -2.40, -2.35, -2.30, -2.25, -2.20, -2.15, -2.10, -2.05, -2.00, -1.95, -1.90, -1.85, -1.80, -1.75, -1.70, -1.65, -1.60, -1.55, -1.50, -1.45, -1.40, -1.35, -1.30, -1.25, -1.20, -1.15, -1.10, -1.05, -1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85, 1.90, 1.95, 2.00, 2.05, 2.10, 2.15, 2.20, 2.25, 2.30, 2.35, 2.40, 2.45, 2.50, 2.55, 2.60, 2.65, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00), ncol = 1)
    scad.goal <- matrix(c(-4.00, -3.95, -3.90, -3.85, -3.80, -3.75, -3.70, -3.65, -3.60, -3.55, -3.50, -3.45, -3.40, -3.35, -3.30, -3.25, -3.20, -3.15, -3.10, -3.05, -3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.40, -2.30, -2.20, -2.10, -2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00), ncol = 1)
    mcp.goal <- matrix(c(-4.000, -3.950, -3.900, -3.850, -3.800, -3.750, -3.700, -3.650, -3.600, -3.550, -3.500, -3.450, -3.400, -3.350, -3.300, -3.250, -3.200, -3.150, -3.100, -3.050, -3.000, -2.925, -2.850, -2.775, -2.700, -2.625, -2.550, -2.475, -2.400, -2.325, -2.250, -2.175, -2.100, -2.025, -1.950, -1.875, -1.800, -1.725, -1.650, -1.575, -1.500, -1.425, -1.350, -1.275, -1.200, -1.125, -1.050, -0.975, -0.900, -0.825, -0.750, -0.675, -0.600, -0.525, -0.450, -0.375, -0.300, -0.225, -0.150, -0.075, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.075, 0.150, 0.225, 0.300, 0.375, 0.450, 0.525, 0.600, 0.675, 0.750, 0.825, 0.900, 0.975, 1.050, 1.125, 1.200, 1.275, 1.350, 1.425, 1.500, 1.575, 1.650, 1.725, 1.800, 1.875, 1.950, 2.025, 2.100, 2.175, 2.250, 2.325, 2.400, 2.475, 2.550, 2.625, 2.700, 2.775, 2.850, 2.925, 3.000, 3.050, 3.100, 3.150, 3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 3.500, 3.550, 3.600, 3.650, 3.700, 3.750, 3.800, 3.850, 3.900, 3.950, 4.000), ncol = 1)

    expect_equal(test_prox_lasso(x, lambda), lasso.goal)
    expect_equal(test_prox_scad(x, lambda, gamma), scad.goal)
    expect_equal(test_prox_mcp(x, lambda, gamma), mcp.goal)
})

test_that("Non-negative prox operators match for non-negative input", {
    lambda <- seq(1, 4)
    x <- matrix(seq(0, 5, 0.05), ncol = 1)

    for (l in lambda) {
        expect_equal(
            test_prox_lasso(x, l),
            test_prox_nnlasso(x, l)
        )

        expect_equal(
            test_prox_scad(x, l),
            test_prox_nnscad(x, l)
        )

        expect_equal(
            test_prox_mcp(x, l),
            test_prox_nnmcp(x, l)
        )
    }
})

test_that("Group lasso proximal operators return correct results", {
    set.seed(32)
    for (grp_prox in c(test_prox_grplasso)) {
        for (rep in 1:10) {
            x <- runif(7)
            # When every element forms a group = lasso
            gp <- as.factor(c(1, 2, 3, 4, 5, 6, 7))
            expect_equal(test_prox_lasso(x, 0), grp_prox(x, gp, 0))
            for (lambda in seq(0, 10, 0.2)) {
                expect_equal(test_prox_lasso(x, lambda), grp_prox(x, gp, lambda))
            }
        }
    }

    for (grp_prox in c(test_prox_grplasso)) {
        # When the all the elements are grouped
        for (rep in 1:10) {
            x <- runif(10)
            gp <- factor(rep(0, 10))
            lambda <- seq(0, 5, 0.2)
            for (l in lambda) {
                scale <- 1 - l / sqrt(sum(x^2))
                if (scale < 0) {
                    scale <- 0
                }
                expect_equal(
                    matrix(scale * x, 10, 1),
                    test_prox_grplasso(x, gp, l)
                )
            }
        }
    }

    # TODO: test for non-negative group lasso

    # When the all the elements are grouped
    for (rep in 1:10) {
        x <- runif(10)
        gp <- factor(rep(0, 10))
        lambda <- seq(0, 1, 0.1)

        for (l in lambda) {
            x_pos <- x * (x > 0)
            scale <- 1 - l / sqrt(sum(x^2))
            expect_equal(
                matrix(scale * x_pos, 10, 1),
                test_prox_nngrplasso(x, gp, l)
            )
            if (scale < 0) {
                scale <- 0
            }
            expect_equal(
                matrix(scale * x_pos, 10, 1),
                test_prox_grplasso(x_pos, gp, l)
            )
        }
    }
})
