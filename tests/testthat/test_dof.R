context("Test degree of freedom")


test_that("DoF of fused lasso", {
    # constant
    x <- c(1, 1, 1, 1)
    expect_equal(1, test_df_orderedfusion(x))

    # fused group at the start
    x <- c(1, 1, 1, 2, 1)
    expect_equal(3, test_df_orderedfusion(x))

    # fused group at the end
    x <- c(1, 2, 1, 1, 1, 1)
    expect_equal(3, test_df_orderedfusion(x))

    # multiple fused groups
    x <- c(1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6)
    expect_equal(6, test_df_orderedfusion(x))

    # no fusion happens
    x <- seq(20)
    expect_equal(20, test_df_orderedfusion(x))
})

test_that("DoF of sparse fused lasso", {
    # constant
    x <- c(1, 1, 1, 1)
    expect_equal(1, test_df_spfusedlasso(x))

    # fused group at the beginning
    x <- c(1, 1, 1, 2, 1)
    expect_equal(3, test_df_spfusedlasso(x))

    # fused group at the end
    x <- c(1, 2, 1, 1, 1, 1)
    expect_equal(3, test_df_spfusedlasso(x))

    # multiple fused groups
    x <- c(1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6)
    expect_equal(6, test_df_spfusedlasso(x))

    # no fusion happens
    x <- seq(20)
    expect_equal(20, test_df_spfusedlasso(x))

    # zeros at the beginning
    x <- c(0, 0, 1, 2, 1, 1)
    expect_equal(3, test_df_spfusedlasso(x))

    # zeros in the middle
    x <- c(1, 0, 0, 2, 1, 1)
    expect_equal(3, test_df_spfusedlasso(x))

    # multiple groups of zeros in the middle
    x <- c(1, 0, 0, 2, 0, 0, 0, 3, 3, 0, 1, 1)
    expect_equal(4, test_df_spfusedlasso(x))

    # zeros in the end
    x <- c(1, 0, 0, 2, 1, 1, 0)
    expect_equal(3, test_df_spfusedlasso(x))
})


test_that("DoF of linear trend filtering", {
    x <- c(1, 1, 1, 1)
    expect_equal(2, test_df_l1gf(x, 1))

    # given any line (knot = 0) it should return DoF = 2
    x <- seq(0, 20, 0.3)
    for (rep in 1:10) {
        alpha <- runif(1)
        beta <- runif(1)
        xx <- alpha * x + beta
        expect_equal(2, test_df_l1gf(x, 1))
    }

    # one knot
    x <- seq(0, 20, 0.3)
    x <- 2 * abs(x - x[37]) - 2
    expect_equal(3, test_df_l1gf(x, 1))

    # three knots
    x <- seq(0, 20, 0.3)
    x <- 2 * abs(x - x[37]) - 2
    x <- 2 * abs(x - x[10]) - 2
    ## plot(x)
    expect_equal(5, test_df_l1gf(x, 1))

    # Find out knots using `diff``
    x <- runif(40)
    change_in_sec_diff <- sum(abs(diff(diff(x))) > 1e-10)
    expect_equal(change_in_sec_diff + 1 + 1, test_df_l1gf(x, 2))
})


test_that("DoF of quadratic trend filtering", {

    # constant
    x <- c(1, 1, 1, 1, 1)
    expect_equal(3, test_df_l1gf(x, 2))

    # given any line it should return DoF = 3
    x <- seq(0, 20, 0.3)
    for (rep in 1:10) {
        alpha <- runif(1)
        beta <- runif(1)
        xx <- alpha * x + beta
        expect_equal(3, test_df_l1gf(x, 2))
    }

    # given any quadratic curve it should return DoF = 3
    x <- seq(0, 20, 0.3)
    for (rep in 1:10) {
        alpha <- runif(1)
        beta <- runif(1)
        c <- runif(1)
        xx <- alpha * x^2 + beta * x + c
        expect_equal(3, test_df_l1gf(x, 2))
    }

    # Find out knots using `diff``
    change_in_sec_diff <- sum(abs(diff(diff(diff(x)))) > 1e-10)

    expect_equal(change_in_sec_diff + 2 + 1, test_df_l1gf(x, 2))
})
