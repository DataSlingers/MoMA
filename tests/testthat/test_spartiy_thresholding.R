context("Thresolding Tests")

test_that("Logging controls work", {
    expect_error(moma_logger_level("BAD LEVEL"))

    moma_logger_level("INFO")
    expect_equal("INFO", moma_logger_level())
    test.points = c(0,0.5, 1, 1.2, 1.5, 2.5,3,3.5)
    x <- seq(-10,10,0.01)
    plot(x,x,type="l")
    lines(x,prox_scad(x,1,3),type="l")
    lines(x,prox_lasso(x,1),type="l",col=1)

    lines(x,prox_mcp(x,1,3),type="l")

    moma_logger_level("MESSAGE")
})
