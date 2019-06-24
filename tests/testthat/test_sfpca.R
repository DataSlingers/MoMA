context("SFPCA tests")
test_that("Equivalent to SVD when no penalty imposed", {
    set.seed(32)

    for (i in 1:30) {
        n <- 17 # set n != p to test bugs
        p <- 23
        X <- matrix(runif(n * p), n)
        sfpca <- sfpca(X)

        svd.result <- svd(X)
        svd.result$u[, 1:4]
        sfpca$u

        expect_equal(norm(svd.result$v[, 1] - sfpca$v), 0)
        expect_equal(norm(svd.result$u[, 1] - sfpca$u), 0)
        expect_equal(svd.result$d[1], sfpca$d[1])
    }
})

test_that("Closed-form solution when Omega = I and no sparsity", {
    set.seed(32)
    n <- 17 # set n != p to test bugs
    p <- 23
    a_u.range <- seq(0, 3, 0.05)
    a_v.range <- seq(0, 3, 0.05)
    for (a_u in a_u.range) {
        for (a_v in a_v.range) {
            for (solver in c("ISTA", "FISTA", "ONESTEPISTA")) {
                # NOTE: We can have one-step ISTA here
                X <- matrix(runif(n * p), n)

                sfpca <- sfpca(X,
                    alpha_u = a_u, alpha_v = a_v, Omega_u = diag(n), Omega_v = diag(p),
                    EPS = 1e-9, MAX_ITER = 1e+5, solver = solver
                )
                svd.result <- svd(X)
                expect_equal(norm(svd.result$v[, 1] - sqrt(1 + a_v) * sfpca$v), 0)
                expect_equal(norm(svd.result$u[, 1] - sqrt(1 + a_u) * sfpca$u), 0)
                expect_equal(svd.result$d[1], sqrt((1 + a_v) * (1 + a_u)) * sfpca$d[1])
            }
        }
    }
})

test_that("Closed-form solution when no sparsity imposed", {
    n <- 17 # set n != p to test bugs
    p <- 23
    set.seed(32)
    X <- matrix(runif(n * p), n)

    # construct p.d. matrix as smoothing matrix
    O_v <- crossprod(matrix(runif(p * p), p, p))
    O_u <- crossprod(matrix(runif(n * n), n, n))

    # set some random alpha's
    # WARNING: running time increases quickly as alpha increases
    a_u.range <- seq(5)
    a_v.range <- seq(5)
    for (a_u in a_u.range) {
        for (a_v in a_v.range) {
            # Cholesky decomposition, note S = I + alpah * Omega
            Lv <- chol(a_v * O_v + diag(p))
            Lu <- chol(a_u * O_u + diag(n))

            svd.result <- svd(t(solve(Lu)) %*% X %*% solve(Lv))
            svd.result.v <- svd.result$v[, 1]
            svd.result.u <- svd.result$u[, 1]

            for (solver in c("ISTA", "FISTA")) {
                # WARNING: One-step ISTA does not pass this test
                res <- sfpca(X,
                    Omega_u = O_u, Omega_v = O_v, alpha_u = a_u, alpha_v = a_v,
                    EPS = 1e-7, MAX_ITER = 1e+5, solver = solver
                )

                # The sfpca solutions and the svd solutions are related by an `L` matrix
                res.v <- Lv %*% res$v
                res.u <- Lu %*% res$u

                # same.direction = 1 if same direction else -1
                same.direction <- ((svd.result$v[, 1][1] * res.v[1]) > 0) * 2 - 1

                # tests
                expect_lte(norm(svd.result$v[, 1] - same.direction * res.v), 1e-5)
                expect_lte(norm(svd.result$u[, 1] - same.direction * res.u), 1e-5)
            }
        }
    }
})

test_that("ISTA and FISTA should yield similar results,
          in the presence of both sparse and smooth penalty", {
    set.seed(332)
    n <- 7 # set n != p to test bugs
    p <- 11
    X <- matrix(runif(n * p), n)

    # generate p.d. matrices
    O_v <- crossprod(matrix(runif(p * p), p, p))
    O_u <- crossprod(matrix(runif(n * n), n, n))

    # run tests
    # NOTE: there's no need to test for large
    # lambda's and alpha's because in those
    # cases u and v are zeros
    cnt <- 0
    for (sp in seq(0, 5, 0.1)) {
        for (sm in seq(0, 5, 0.1)) {
            # TODO: Add "L1TRENDFILTERING"
            for (sptype in c("LASSO", "SCAD", "MCP", "ORDEREDFUSED")) {
                ista <- sfpca(X,
                    Omega_u = O_u, Omega_v = O_v, alpha_u = sp, alpha_v = sp,
                    lambda_u = sm, lambda_v = sm, P_u = "LASSO", P_v = sptype,
                    EPS = 1e-14, MAX_ITER = 1e+3, solver = "ISTA", EPS_inner = 1e-9
                )
                fista <- sfpca(X,
                    Omega_u = O_u, Omega_v = O_v, alpha_u = sp, alpha_v = sp,
                    lambda_u = sm, lambda_v = sm, P_u = "LASSO", P_v = sptype,
                    EPS = 1e-6, MAX_ITER = 1e+3, solver = "FISTA", EPS_inner = 1e-9
                )

                # WARNING: We observe if zero appears in either v or u, ista and fista
                # might not give identical results.
                # Maybe they will both eventually go to the same point, but ista slows
                # down a lot before it reaches it and consequently meets the stopping criterion.
                if (sum(ista$v[, 1] == 0.0) == 0
                && sum(fista$v[, 1] == 0.0) == 0) {
                    expect_lte(sum((ista$v[, 1] - fista$v[, 1])^2), 1e-6)
                    expect_lte(sum((ista$u[, 1] - fista$u[, 1])^2), 1e-6)
                }
            }
        }
    }
})
