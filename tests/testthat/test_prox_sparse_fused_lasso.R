context("Sparse fused lasso tests")

test_that("Same results as the `flsa` package", {
    set.seed(43)
    if (requireNamespace("flsa")) {
        library(flsa)
        pset <- seq(2, 8)
        for (p in pset) {
            for (i in 1:20) {
                x <- 10 * runif(p)
                for (lambda in seq(0, 2, 0.1)) {
                    for (lambda2 in seq(0, 2, 0.2)) {
                        expect_equal(
                            test_prox_spfusedlasso(x, lambda, lambda2 = lambda2),
                            matrix(flsaGetSolution(flsa(x), lambda2 = lambda, lambda1 = lambda2))
                        )
                    }
                }
            }
        }
    }
})
