context("ISTA Tests")
test_that("Equivalent to SVD when no penalty imposed", {
    set.seed(32)
    n <- 10 # set n != p to test bugs
    p <- 6

    set.seed(32)
    X = matrix(runif(n*p),n)

    for(P_u in list("lasso","scad","mcp")){
        for(P_v in list("lasso","scad","mcp")){
            sfpca <- sfpca(X,
                           alpha_u=0,alpha_v=0,
                           lambda_u=0,lambda_v=0,P_u=P_u,P_v=P_v,
                           EPS=1e-9,MAX_ITER = 1e+5)
            svd.result <- svd(X)
            svd.result$u[,1:4]
            sfpca$u
            expect_equal(norm(svd.result$v[,1]-sfpca$v),0)
            expect_equal(norm(svd.result$u[,1]-sfpca$u),0)
            expect_equal(svd.result$d[1],sfpca$d[1]);
        }
    }
})

test_that("Equivalent to SVD when Omega = I and no sparsity",{
    set.seed(32)
    n <- 10 # set n != p to test bugs
    p <- 6
    a_u.range <- seq(0,3,0.3)
    a_v.range <- seq(0,3,0.4)
    for(a_u in a_u.range){
        for(a_v in a_v.range){
            for(solver in c("ISTA","FISTA","ONESTEPISTA")){
                # NOTE: We can have one-step ISTA here
                a_v <- 1
                set.seed(32)
                X = matrix(runif(n*p),n)

                sfpca <- sfpca(X,
                               alpha_u=a_u,alpha_v=a_v,Omega_u=diag(n),Omega_v=diag(p),
                               lambda_u=0,lambda_v=0,P_u="LASSO",P_v="LASSO",
                               EPS=1e-9,MAX_ITER = 1e+5,solver = solver)
                svd.result <- svd(X)
                expect_equal(norm(svd.result$v[,1] - sqrt(1 + a_v) * sfpca$v),0)
                expect_equal(norm(svd.result$u[,1] - sqrt(1 + a_u) * sfpca$u),0)
                expect_equal(svd.result$d[1],sqrt((1 + a_v) * (1 + a_u)) * sfpca$d[1]);
            }

        }
    }
})

test_that("Closed form solution when no sparsity imposed",{
    n <- 30 # set n != p to test bugs
    p <- 40
    set.seed(32)
    X = matrix(runif(n*p),n)

    # construct p.d. matrix as smoothing matrix
    mat1 <- matrix(runif(p*p),p,p)
    mat2 <- matrix(runif(n*n),n,n)
    O_v = t(mat1) %*% mat1
    O_u = t(mat2) %*% mat2

    # set some random alpha's
    # WARNING: running time increases quickly as alpha increases
    a_u.range <- c(1,2,3,4)
    a_v.range <- c(5)
    for(a_u in a_u.range){
        for(a_v in a_v.range){
            for(solver in c("ISTA","FISTA")){
                # WARNING: One-step ISTA does not pass this test

                # Cholesky decomposition, note S = I + alpah * Omega
                Lv = chol(a_v * O_v + diag(p))
                Lu = chol(a_u * O_u + diag(n))

                svd.result <- svd(t(solve(Lu)) %*% X %*% solve(Lv))
                svd.result.v = svd.result$v[,1]
                svd.result.u = svd.result$u[,1]

                res <- sfpca(X,
                              Omega_u=O_u,Omega_v=O_v,alpha_u=a_u,alpha_v=a_v,
                              EPS=1e-9,MAX_ITER=1e+5,solve=solver)
                # The sfpca solutions and the svd solutions are related by an `L` matrix
                res.v = Lv %*% res$v
                res.u = Lu %*% res$u

                # same.direction = 1 if same direction else -1
                same.direction = ((svd.result$v[,1][1] * res.v[1]) > 0) * 2 -1

                # tests
                expect_lte(norm(svd.result$v[,1] - same.direction * res.v),1e-4)
                expect_lte(norm(svd.result$u[,1] - same.direction * res.u),1e-5)
            }
        }
    }
})

test_that("ISTA and FISTA should yield similar results",{
    set.seed(302)
    n <- 6 # set n != p to test bugs
    p <- 10
    X = matrix(runif(n*p),n)

    # generate p.d. matrix
    mat1 <- matrix(runif(p*p),p,p)
    mat2 <- matrix(runif(n*n),n,n)
    O_v = t(mat1) %*% mat1
    O_u = t(mat2) %*% mat2

    # run algorithms
    for(sp in c(0,1,2,3,4)){
        for(sptype in c("LASSO","SCAD","MCP")){
            for(sm in seq(0,0.3,0.15)){
                ista <- sfpca(X,
                              Omega_u=O_u,Omega_v=O_v,alpha_u=sp,alpha_v=sp,
                              lambda_u=sm,lambda_v=sm,P_u="LASSO",P_v=sptype,
                              EPS=1e-6,MAX_ITER = 1e+6,solve="ISTA")
                fista <- sfpca(X,
                               Omega_u=O_u,Omega_v=O_v,alpha_u=sp,alpha_v=sp,
                               lambda_u=sm,lambda_v=sm,P_u="LASSO",P_v=sptype,
                               EPS=1e-6,MAX_ITER = 1e+5,solve="FISTA")
                expect_lte(sum((ista$v[,1]-fista$v[,1])^2),1e-6)
                expect_lte(sum((ista$u[,1]-fista$u[,1])^2),1e-6)
            }
        }
    }
})
