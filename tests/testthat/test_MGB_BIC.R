context("`cpp_multirank_BIC_grid_search` as greedy BIC")

set.seed(123)
X <- matrix(runif(12), 3, 4) * 10

# They are mutually prime. This is useful for testing
# the size of the 5D list
n_alpha_u <- 7
n_alpha_v <- 5
n_lambda_u <- 3
n_lambda_v <- 2
rank <- 3

# Penalties should not be too high. Otherwise
# u and v becomes vectors of zeros, and defaltion
# is thus trivial.
alpha_u <- seq(0, 3, length.out = n_alpha_u)
alpha_v <- seq(0, 3, length.out = n_alpha_v)
lambda_u <- seq(0, 3, length.out = n_lambda_u)
lambda_v <- seq(0, 3, length.out = n_lambda_v)

# public_arglist_wo_rank_and_selection does not specify two things: rank
# and selection strategy for each parameter. They
# are specified in each test case.
public_arglist_wo_rank_and_selection <- c(
    list(
        X = X,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(3), Omega_v = second_diff_mat(4),
        lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = add_default_prox_args(lasso()),
        prox_arg_list_v = add_default_prox_args(lasso())
    ),
    moma_pg_settings()
)

test_that("Returns correct vectors for chosen parameters", {
    result <- do.call(
        cpp_multirank_BIC_grid_search,
        c(
            public_arglist_wo_rank_and_selection,
            list(
                select_scheme_alpha_u = 0, # grid
                select_scheme_alpha_v = 0, # grid
                select_scheme_lambda_u = 1, # bic
                select_scheme_lambda_v = 1 # bic
            )
        )
    )

    # Loop order in C++ is (outmost) au, lu, av, lv (innermost)
    expect_equal(dim(result), c(7, 1, 5, 1, 1))

    for (i in 1:n_alpha_u) {
        for (j in 1:n_alpha_v) {
            res_i_j <- get_5Dlist_elem(result, i, 1, j, 1)[[1]]
            opt_alpha_u <- res_i_j$u$alpha
            opt_alpha_v <- res_i_j$v$alpha
            opt_lambda_u <- res_i_j$u$lambda
            opt_lambda_v <- res_i_j$v$lambda

            arglist_unary_penalty <- modifyList(
                public_arglist_wo_rank_and_selection,
                list(
                    alpha_u = opt_alpha_u, alpha_v = opt_alpha_v,
                    lambda_u = opt_lambda_u, lambda_v = opt_lambda_v
                )
            )

            result_directcall <- do.call(
                cpp_moma_multi_rank,
                arglist_unary_penalty
            )

            expect_equal(result_directcall$u, res_i_j$u$vector)
            expect_equal(result_directcall$v, res_i_j$v$vector)
        }
    }
})
