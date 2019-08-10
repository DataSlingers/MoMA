context("BIC tests")

bic_lasso <- function(y, y_est) {
    p <- length(y)
    res <- norm(as.matrix(y - y_est), "2")
    df <- sum(y_est != 0)
    bic <- log(res * res / p) + log(p) / p * df
    return(bic)
}

test_that("Test for lasso BIC", {
    y <- c(1, 2, 3)
    y_est <- c(2, 2, 2)
    p <- length(y)
    expect_equal(test_BIC(
        y, y_est,
        "ISTA",
        0, second_diff_mat(p),
        0, add_default_prox_args(lasso()),
        p
    ), bic_lasso(y, y_est))
})
