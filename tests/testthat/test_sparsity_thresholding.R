context("Thresolding Tests")

test_that("Input validity test", {
    # TODO
})
test_that("When lambda=0, should return exactly the input", {

    lambda <- 0
    gamma <- 3
    x <- seq(-4,4,0.01)
    # run those operators
    scadx <- prox_scad(x,lambda,3)
    mcpx <- prox_mcp(x,lambda,3)
    lassox <- prox_lasso(x,lambda)

    # tests
    expect_equal(norm(scadx - x),0)
    expect_equal(norm(mcpx - x),0)
    expect_equal(norm(lassox - x),0)

})

test_that("When lambda=3", {
    # TODO
})

test_that("Group lasso proximal operator test", {
    # TODO
})

