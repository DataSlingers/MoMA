context("ISTA Tests")
test_that("Equivalent to SVD when no penalty imposed", {
    n <- 10 # set n != p to test bugs
    p <- 6
    O_v <- diag(p)
    O_u <- diag(n)
    set.seed(32)
    X = matrix(runif(n*p),n)


    sfpca <- sfpca(X,
                   O_u,O_v, 0,0,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5)
    svd.result <- svd(X)
    expect_equal(norm(svd.result$v[,1]-sfpca$v),0)
    expect_equal(norm(svd.result$u[,1]-sfpca$u),0)
    expect_equal(svd.result$d[1],sfpca$d);
})

test_that("Equivalent to SVD when Omega = I and no sparsity",{
    n <- 10 # set n != p to test bugs
    p <- 6
    O_v <- diag(p)
    O_u <- diag(n)
    set.seed(32)
    X = matrix(runif(n*p),n)

    sfpca <- sfpca(X,
                   O_u,O_v,1/n,1/p,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5)
    svd.result <- svd(X)
    expect_equal(norm(svd.result$v[,1]-sfpca$v),0)
    expect_equal(norm(svd.result$u[,1]-sfpca$u),0)
    expect_equal(svd.result$d[1],sfpca$d);
})


test_that("Closed form solution when no sparsity imposed",{
    n <- 30 # set n != p to test bugs
    p <- 40
    set.seed(32)
    X = matrix(runif(n*p),n)

    # construct p.d. matrix as smoothing matrix
    Xv <- matrix(runif(p*p),p,p)
    Xu <- matrix(runif(n*n),n,n)
    O_v = t(Xv) %*% Xv
    O_u = t(Xu) %*% Xu

    # Cholesky decomposition, note S=I+alpah*Omega, alpha set to 1 here
    Lv = chol(O_v+diag(p))
    Lu = chol(O_u+diag(n))

    svd.result <- svd(t(solve(Lu)) %*% X %*% solve(Lv))
    svd.result.v = svd.result$v[,1] / sqrt(sum(svd.result$v[,1]^2))
    svd.result.u = svd.result$u[,1] / sqrt(sum(svd.result$u[,1]^2))
    ista <- sfpca(X,
                   O_u,O_v,1,1,
                   EPS=1e-9,MAX_ITER = 1e+5,solve="ISTA")
    fista <- sfpca(X,
                   O_u,O_v,1,1,
                   EPS=1e-9,MAX_ITER = 1e+5,solve="FISTA")

    # know that sfpca always return norm 1 vector
    ista.v = Lv %*% ista$v / norm(Lv %*% ista$v,"E")
    ista.u = Lu %*% ista$u / norm(Lu %*% ista$u,"E")
    fista.v = Lv %*% fista$v / norm(Lv %*% fista$v,"E")
    fista.u = Lu %*% fista$u / norm(Lu %*% fista$u,"E")

    # different up to sign
    expect_lte(norm(svd.result$v[,1] - ista.v),1e-5)
    expect_lte(norm(svd.result$u[,1] - ista.u),1e-5)
    expect_lte(norm(svd.result$v[,1] - fista.v),1e-5)
    expect_lte(norm(svd.result$u[,1] - fista.u),1e-6)
})

test_that("ISTA and FISTA should yield similar results",{
    n <- 6 # set n != p to test bugs
    p <- 10
    X = matrix(runif(n*p),n)

    # generate p.d. matrix
    Xv <- matrix(runif(p*p),p,p)
    Xu <- matrix(runif(n*n),n,n)
    O_v = t(Xv) %*% Xv
    O_u = t(Xu) %*% Xu

    # run algorithms
    svd.result <- svd(t(solve(O_u)) %*% X %*% solve(O_v))
    ista <- sfpca(X,
                   O_u,O_v,1,1,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5,solve="ISTA")
    fista <- sfpca(X,
                   O_u,O_v,1,1,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5,solve="FISTA")

    # tests
    expect_equal(ista$v[,1],fista$v[,1])
})
