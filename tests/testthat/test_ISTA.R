
#-------------------
# test_that
#-------------------
test_that("Input validity test",{
})

test_that("Equal to SVD when no penalty", {
    sfpca <- sfpca(X,
                   O_u,O_v, 0,0,
                   lambda_u=0,lambda_v=0,"LASSO","LASSO",
                   gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5,)
    svdd <- svd(X)
    expect_equal(sum((svdd$v[,1]-sfpca$v)^2),0)
    expect_equal(sum((svdd$u[,1]-sfpca$u)^2),0)
    expect_equal(svdd$d[1],sfpca$d);
    expect_error(moma_logger_level("BAD LEVEL"))
})

test_that("Equal to analytic solution when only roughness penalty n", {
    # TODO
})

test_that("Closed form solution when no sparsity",{
    # TODO
})
