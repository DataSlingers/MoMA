context("Test BIC-grid-mixed search")


test_that("BIC search returns correct-sized grid", {
    set.seed(123)
    X <- matrix(runif(12),3,4)

    n_au=7; n_av=5; n_lu=3; n_lv=2; # mutually prime
    alpha_u=seq(0.3,1,length.out = n_au)
    alpha_v=seq(0.3,1,length.out = n_av)
    lambda_u=seq(0.3,1,length.out = n_lu)
    lambda_v=seq(0.3,1,length.out = n_lv)


    arg_list <- list(
        X,
        alpha_u=alpha_u, alpha_v=alpha_v,
        Omega_u=second_diff_mat(3), Omega_v=second_diff_mat(4),
        lambda_u=lambda_u, lambda_v=lambda_v,
        prox_arg_list_u=add_default_prox_args(lasso()), prox_arg_list_v=add_default_prox_args(empty()),
        EPS=1e-6, MAX_ITER=1e+4, EPS_inner=1e-6, MAX_ITER_inner=1e+4, solver="ISTA"
    )

    # Case 1: four grid requests
    result <- do.call(testnestedBIC,
                      c(arg_list,
                        list(bicau=0,
                             bicav=0,
                             biclu=0,
                             biclv=0)))

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv = sapply(result, function(x) x$v$lambda)
    av = sapply(result, function(x) x$v$alpha)
    lu = sapply(result, function(x) x$u$lambda)
    au = sapply(result, function(x) x$u$alpha)

    expect_equal(lv, rep(lambda_v, n_au*n_lu*n_av, each=1))
    expect_equal(av, rep(alpha_v,  n_au*n_lu,      each=n_lv))
    expect_equal(lu, rep(lambda_u, n_au,           each=n_lv*n_av))
    expect_equal(au, rep(alpha_u,  1,              each=n_lv*n_av*n_lu))


    # Case 2: three grid requests, one BIC requet
    result2 <- do.call(testnestedBIC,
                      c(arg_list,
                        list(bicau=0,
                             bicav=0,
                             biclu=0,
                             biclv=1)))

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    # lv = sapply(result2, function(x) x$v$lambda)
    av = sapply(result2, function(x) x$v$alpha)
    lu = sapply(result2, function(x) x$u$lambda)
    au = sapply(result2, function(x) x$u$alpha)

    # expect_equal(lv, rep(lambda_v, n_au*n_lu*n_av, each=1))
    expect_equal(av, rep(alpha_v,  n_au*n_lu, each=1))
    expect_equal(lu, rep(lambda_u, n_au,      each=n_av))
    expect_equal(au, rep(alpha_u,  1,         each=n_lu*n_av))
})
