context("`cpp_multirank_BIC_grid_search` as grid-and-BIC")

set.seed(123)
X <- matrix(runif(12), 3, 4) * 10

# They are mutually prime. This is useful for testing
# the size of the 5D list
n_alpha_u <- 7
n_alpha_v <- 5
n_lambda_u <- 3
n_lambda_v <- 2

# Penalties should not be too high. Otherwise
# u and v becomes vectors of zeros, and defaltion
# is thus trivial.
alpha_u <- seq(0, 1, length.out = n_alpha_u)
alpha_v <- seq(0, 2, length.out = n_alpha_v)
lambda_u <- seq(0, 3, length.out = n_lambda_u)
lambda_v <- seq(0, 4, length.out = n_lambda_v)

# public_arg_list_wo_rank_and_selection does not specify two things: rank
# and selection strategy for each parameter. They
# are specified in each test case.
public_arg_list_wo_rank_and_selection <- c(
    list(
        X = X,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(3), Omega_v = second_diff_mat(4),
        lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = add_default_prox_args(lasso()), prox_arg_list_v = add_default_prox_args(lasso())
    ),
    moma_pg_settings()
)

# Tests for multirank + BIC
test_that("Returns correct-sized grid: four grid requests", {
    # Case 1: four grid requests
    result <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 0, # grid
                select_scheme_alpha_v = 0, # grid
                select_scheme_lambda_u = 0, # grid
                select_scheme_lambda_v = 0 # grid
            )
        )
    )

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result, function(x) x$v$lambda)
    av <- sapply(result, function(x) x$v$alpha)
    lu <- sapply(result, function(x) x$u$lambda)
    au <- sapply(result, function(x) x$u$alpha)

    expect_equal(dim(result), c(n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v, 1))

    expect_equal(lv, rep(lambda_v, n_alpha_u * n_lambda_u * n_alpha_v, each = 1))
    expect_equal(av, rep(alpha_v, n_alpha_u * n_lambda_u, each = n_lambda_v))
    expect_equal(lu, rep(lambda_u, n_alpha_u, each = n_lambda_v * n_alpha_v))
    expect_equal(au, rep(alpha_u, 1, each = n_lambda_v * n_alpha_v * n_lambda_u))
})

test_that("Returns correct-sized grid: three grid requests", {
    # BIC on lambda_v
    result2 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 0, # grid
                select_scheme_alpha_v = 0, # grid
                select_scheme_lambda_u = 0, # grid
                select_scheme_lambda_v = 1
            )
        )
    )

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    av <- sapply(result2, function(x) x$v$alpha)
    lu <- sapply(result2, function(x) x$u$lambda)
    au <- sapply(result2, function(x) x$u$alpha)

    expect_equal(dim(result2), c(n_alpha_u, n_lambda_u, n_alpha_v, 1, 1))

    expect_equal(av, rep(alpha_v, n_alpha_u * n_lambda_u, each = 1))
    expect_equal(lu, rep(lambda_u, n_alpha_u, each = n_alpha_v))
    expect_equal(au, rep(alpha_u, 1, each = n_lambda_u * n_alpha_v))

    # BIC on alpha_u
    result2 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 1,
                select_scheme_alpha_v = 0, # grid
                select_scheme_lambda_u = 0, # grid
                select_scheme_lambda_v = 0
            )
        )
    ) # grid

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result2, function(x) x$v$lambda)
    av <- sapply(result2, function(x) x$v$alpha)
    lu <- sapply(result2, function(x) x$u$lambda)

    expect_equal(dim(result2), c(1, n_lambda_u, n_alpha_v, n_lambda_v, 1))

    expect_equal(lv, rep(lambda_v, n_lambda_u * n_alpha_v, each = 1))
    expect_equal(av, rep(alpha_v, n_lambda_u, each = n_lambda_v))
    expect_equal(lu, rep(lambda_u, 1, each = n_alpha_v * n_lambda_v))
})

test_that("Returns correct-sized grid: two grid requests on u", {
    # Case 3: two grid requests, both on u side, and two BIC
    result3 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 0, # grid
                select_scheme_lambda_u = 0, # grid
                select_scheme_alpha_v = 1,
                select_scheme_lambda_v = 1
            )
        )
    )

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result3, function(x) x$v$lambda)
    av <- sapply(result3, function(x) x$v$alpha)
    lu <- sapply(result3, function(x) x$u$lambda)
    au <- sapply(result3, function(x) x$u$alpha)

    expect_equal(dim(result3), c(n_alpha_u, n_lambda_u, 1, 1, 1))

    expect_equal(lu, rep(lambda_u, n_alpha_u, each = 1))
    expect_equal(au, rep(alpha_u, 1, each = n_lambda_u))
})

test_that("Returns correct-sized grid: two grid requests on different sides", {
    # Case 4: two grid requests, both on u side, and two BIC
    result4 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 1,
                select_scheme_lambda_u = 0, # grid
                select_scheme_alpha_v = 1,
                select_scheme_lambda_v = 0
            )
        )
    ) # grid

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result4, function(x) x$v$lambda)
    av <- sapply(result4, function(x) x$v$alpha)
    lu <- sapply(result4, function(x) x$u$lambda)
    au <- sapply(result4, function(x) x$u$alpha)

    expect_equal(dim(result4), c(1, n_lambda_u, 1, n_lambda_v, 1))

    expect_equal(lv, rep(lambda_v, n_lambda_u, each = 1))
    expect_equal(lu, rep(lambda_u, 1, each = n_lambda_v))
})

test_that("Returns correct-sized grid: one grid", {
    # Case 5: one grid requests, both on u side, and two BIC
    result4 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 1,
                select_scheme_lambda_u = 1,
                select_scheme_alpha_v = 1,
                select_scheme_lambda_v = 0
            )
        )
    ) # grid

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result4, function(x) x$v$lambda)
    av <- sapply(result4, function(x) x$v$alpha)
    lu <- sapply(result4, function(x) x$u$lambda)
    au <- sapply(result4, function(x) x$u$alpha)

    expect_equal(dim(result4), c(1, 1, 1, n_lambda_v, 1))

    expect_equal(lv, rep(lv, 1, each = 1))
})

test_that("Returns correct-sized grid: four BIC search", {
    # Case 6: one grid requests, both on u side, and two BIC
    result4 <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arg_list_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 1,
                select_scheme_lambda_u = 1,
                select_scheme_alpha_v = 1,
                select_scheme_lambda_v = 1
            )
        )
    )

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    lv <- sapply(result4, function(x) x$v$lambda)
    av <- sapply(result4, function(x) x$v$alpha)
    lu <- sapply(result4, function(x) x$u$lambda)
    au <- sapply(result4, function(x) x$u$alpha)

    expect_equal(dim(result4), c(1, 1, 1, 1, 1))

    expect_equal(lv, rep(lv, 1, each = 1))
})

test_that("`cpp_multirank_BIC_grid_search` receives a vector of length 0", {
    arglist_w_empty_penalty <- c(
        modifyList(
            public_arg_list_wo_rank_and_selection,
            list(lambda_u = vector())
        ),
        list(
            select_scheme_alpha_u = 1,
            select_scheme_lambda_u = 1,
            select_scheme_alpha_v = 1,
            select_scheme_lambda_v = 1
        )
    )

    expect_error(
        do.call(cpp_multirank_BIC_grid_search, arglist_w_empty_penalty),
        "Please specify all four parameters"
    )
})
