context("ISTA Tests")
#-------------------
# Generate data
#-------------------
n <- 10 # set n != p to test bugs
p <- 6
O_v <- diag(p)
O_u <- diag(n)
set.seed(32)
X = matrix(runif(n*p),n)


#-------------------
# test_that
#-------------------

test_that("Equivalent to SVD when no penalty imposed", {

    sfpca <- sfpca(X,
                   O_u,O_v, 0,0,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5)
    svdd <- svd(X)
    expect_equal(norm(svdd$v[,1]-sfpca$v),0)
    expect_equal(norm(svdd$u[,1]-sfpca$u),0)
    expect_equal(svdd$d[1],sfpca$d);
})

test_that("Equivalent to SVD when Omega = I and no sparsity",{
    sfpca <- sfpca(X,
                   O_u,O_v,1/n,1/p,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5)
    svdd <- svd(X)
    expect_equal(norm(svdd$v[,1]-sfpca$v),0)
    expect_equal(norm(svdd$u[,1]-sfpca$u),0)
    expect_equal(svdd$d[1],sfpca$d);
})

test_that("Closed form solution when no sparsity imposed",{
    # Xv <- matrix(runif(p*p),p,p)
    # Xu <- matrix(runif(n*n),n,n)

    # O_v = t(Xv) %*% Xv
    # O_u = t(Xu) %*% Xu
    # Lv = chol(O_v)
    # Lu = chol(O_u)
    # svdd <- svd(t(solve(Lu)) %*% X %*% solve(Lv))
    # sfpca <- sfpca(X,
    #                O_u,O_v,1/n,1/p,
    #                lambda_u=0,lambda_v=0,"LASSO","LASSO",
    #                gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5)
    # expect_equal(norm(svdd$v[,1]-Lv %*% sfpca$v),0)
    # expect_equal(norm(svdd$u[,1]-Lu %*% sfpca$u),0)
    # expect_equal(svdd$d[1],sfpca$d);
    # TODO
})
