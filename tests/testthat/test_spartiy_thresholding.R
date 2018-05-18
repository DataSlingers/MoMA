context("Thresolding Tests")

test_that("Gamma validation", {
    x <- c(1,2,3)
    expect_error(prox_scad(x,1,1.999))
    expect_error(prox_mcp(x,1,0.999))

})

test_that("Correct calculation", {
    x <- seq(-5,5,0.3)
    expect_equal(prox_scad(x,1,3),t(t(as.vector(c(-5.0, -4.7, -4.4, -4.1, -3.8, -3.5, -3.2, -2.8, -2.2,
                                                   -1.6, -1.0, -0.7, -0.4, -0.1,  0.0,  0.0,  0.0,  0.0,
                                                   0.0,  0.0,  0.0,  0.3,  0.6,  0.9,  1.4,  2.0,  2.6,
                                                   3.1,  3.4,  3.7,  4.0,  4.3,  4.6,  4.9)))))
    expect_equal(prox_mcp(x,1,3),t(t(as.vector(-5.00, -4.70, -4.40, -4.10, -3.80, -3.50, -3.20,
                                               -2.85, -2.40, -1.95, -1.50, -1.05, -0.60, -0.15,
                                               0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                                               0.45,  0.90,  1.35,  1.80,  2.25,  2.70,  3.10,
                                               3.40,  3.70,  4.00,  4.30,  4.60,  4.90))))

    x <- seq(-4,4,0.01)
    plot(x,x,type="l")
    lines(x,prox_scad(x,1,3),type="l")

    lines(x,prox_mcp(x,1,3),type="l")
    lines(x,prox_lasso(x,1),type="l",col=1)
    })
