context("Unordered Fusion Lasso Tests")

# Special case 1:
# Large lambda: lambdas in each group are
# the means of everything in the group
test_that("Find means of everything when lambda is large enough and the graph is fully connected", {
    set.seed(43)
    rep <- 10
    large.lambda <- 10000
    for (p in c(5, 30, 40, 100)) {
        y <- runif(p)

        # A fully connected graph with random weights
        w <- matrix(rep(0, p * p), p, byrow = T)
        for (i in 1:p - 1) {
            for (j in (i + 1):p) {
                w[i, j] <- runif(1) + 1
            }
        }

        for (i in 1:rep) {
            res.AMA <- test_prox_fusion(y, large.lambda, w, ADMM = FALSE, acc = FALSE)
            res.AMA.acc <- test_prox_fusion(y, large.lambda, w, ADMM = FALSE, acc = TRUE)
            res.ADMM <- test_prox_fusion(y, large.lambda, w, ADMM = TRUE, acc = FALSE)
            for (j in c(
                res.AMA,
                res.AMA.acc,
                res.ADMM
            )) {
                expect_equal(j, mean(y))
            }
        }
    }
})

test_that("Find means of connected components when lambda is large enough", {
    set.seed(33)
    rep <- 10
    large.lambda <- 1000
    for (p in c(9, 99)) {

        # A weight matrix of 3 connected components
        w <- matrix(rep(0, p * p), p, byrow = T)
        num.comp.nodes <- p / 3 # the number of nodes in each component
        for (i in seq(1, p, num.comp.nodes)) {
            w[i, (i + 1):(i + num.comp.nodes - 1)] <- 1 + runif(num.comp.nodes - 1)
        }

        for (j in 1:rep) {
            y <- 10 * runif(p)
            res.AMA.unacc <- test_prox_fusion(y, large.lambda, w, ADMM = FALSE, acc = TRUE)
            res.AMA.acc <- test_prox_fusion(y, large.lambda, w, ADMM = FALSE, acc = TRUE)
            res.ADMM <- test_prox_fusion(y, large.lambda, w, ADMM = TRUE, acc = FALSE)
            for (res in list(
                res.AMA.unacc,
                res.AMA.acc,
                res.ADMM
            )) {
                for (i in 1:p) {
                    start <- ((i - 1) %/% num.comp.nodes) * num.comp.nodes + 1
                    end <- start + num.comp.nodes - 1
                    expect_lte(abs(res[i] - mean(y[start:end])), 1e-6)
                }
            }
        }
    }
})

# Special case 2:
# For every lambda: when only w_ij = 1 (j = i+1)
# it becomes ordered fused lasso
test_that("Ordered fused lasso when w_ij = 1 all j = i+1", {
    set.seed(44)
    rep <- 20
    for (p in c(3, 20, 100)) {

        # A chained graph
        w <- matrix(rep(0, p * p), p, byrow = T)
        for (i in 1:(p - 1)) {
            w[i, i + 1] <- 1
        }

        for (i in 1:rep) {
            y <- 10 * runif(p)
            err.AMA.unacc <- norm(test_prox_fusedlassodp(y, 1) - test_prox_fusion(y, 1, w, ADMM = FALSE, acc = TRUE))
            err.AMA.acc <- norm(test_prox_fusedlassodp(y, 1) - test_prox_fusion(y, 1, w, ADMM = FALSE, acc = TRUE))
            err.ADMM <- norm(test_prox_fusedlassodp(y, 1) - test_prox_fusion(y, 1, w, ADMM = TRUE, acc = FALSE))
            for (err in c(
                err.AMA.unacc,
                err.AMA.acc,
                err.ADMM
            )) {
                expect_lte(err, 1e-4)
            }
        }
    }
})

# Special case 3:
# For evry lambda: when all i,j w_ij = 1
# the solution path contains no split.

# Find the weight vector required by the
# `cvxclustr` package
mat_to_vec <- function(my.w, p) {
    cnt <- 1
    cvx.w <- vector(mode = "numeric", length = p * (p - 1) / 2)
    for (jj in 1:(p - 1)) {
        for (ii in (jj + 1):p) {
            cvx.w[cnt] <- my.w[jj, ii]
            cnt <- cnt + 1
        }
    }
    return(cvx.w)
}

test_that("Unweighted and fully connected graph, i.e., w_ij = 1 for all i, j", {
    if (requireNamespace("cvxclustr")) {
        set.seed(33)
        rep <- 20

        for (p in c(30, 100)) {

            # A weight matrix where w_ij = 1, all i, j
            w <- matrix(rep(1, p * p), p, byrow = T)
            for (i in 1:(p - 1)) {
                for (j in (i + 1):p) {
                    w[i, j] <- runif(1) + 1
                }
            }

            for (i in 1:rep) {
                y <- 10 * runif(p)
                y <- t(matrix(y))

                # The cvxclustr package stores weights as a vector.
                # For Gaussian kernel, wij = exp(-phi ||X[,i]-X[,j]||^2)
                # So phi = 0 makes a fully connected and unweighted graph
                for (lambda in seq(0, 1, 0.13)) {
                    cvx.result <- cvxclustr::cvxclust(y, mat_to_vec(w, p), lambda, method = "admm", tol = 1e-10)$U[[1]]

                    admm <- t(matrix(test_prox_fusion(y, lambda, w, ADMM = TRUE, acc = FALSE)))
                    ama <- t(matrix(test_prox_fusion(y, lambda, w, ADMM = FALSE, acc = FALSE)))
                    ama.acc <- t(matrix(test_prox_fusion(y, lambda, w, ADMM = FALSE, acc = TRUE)))
                }

                err.AMA.unacc <- norm(cvx.result - ama)
                err.AMA.acc <- norm(cvx.result - ama.acc)
                err.ADMM <- norm(cvx.result - admm)
                for (err in c(
                    err.AMA.unacc,
                    err.AMA.acc,
                    err.ADMM
                )) {
                    expect_lte(err, 1e-6)
                }
            }
        }
    }
})
