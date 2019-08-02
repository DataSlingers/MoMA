context("SLOPE")

proxSortedL1 <- function(x, l) {

    # Modification: Choose HB type lambda
    lambda <- vector(mode = "numeric", length = length(x))

    for (i in 1:length(x)) {
        lambda[i] <- qnorm(1 - i * 0.05 / 2 / length(x))
    }

    lambda <- lambda * l

    # Compare to the SLOPE package on CRAN
    result <- SLOPE::prox_sorted_L1(x, lambda, method = "c")

    return(result)
}

test_that("Compared to the SLOPE package", {
    if (requireNamespace("SLOPE")) {
        set.seed(123)
        reps <- 100
        for (i in 1:reps) {
            x <- runif(10)
            for (lambda in seq(0, 3, 0.2)) {
                expect_equal(
                    test_prox_slope(x, lambda),
                    as.matrix(proxSortedL1(x, lambda))
                )
            }
        }
    }
})
