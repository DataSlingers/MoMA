context("`cpp_multirank_BIC_grid_search` as solving a single problem")

set.seed(123)
X <- matrix(runif(12), 3, 4) * 10

# They are mutually prime. This is useful for testing
# the size of the 5D list
n_alpha_u <- 7
n_alpha_v <- 5
n_lambda_u <- 5
n_lambda_v <- 5
rank <- 1

# Penalties should not be too high. Otherwise
# u and v becomes vectors of zeros, and defaltion
# is thus trivial.
alpha_u <- seq(0, 10, length.out = n_alpha_u)
alpha_v <- seq(0, 10, length.out = n_alpha_v)
lambda_u <- seq(0, 10, length.out = n_lambda_u)
lambda_v <- seq(0, 10, length.out = n_lambda_v)

# public_arglist_wo_penalty does not specify
public_arglist_wo_penalty <- c(
    list(
        X = X,
        #  alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(3), Omega_v = second_diff_mat(4),
        #  lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = add_default_prox_args(lasso()),
        prox_arg_list_v = add_default_prox_args(lasso()),
        rank = 1
    ),
    moma_pg_settings()
)

# Tests for greedy BIC
test_that("`cpp_multirank_BIC_grid_search` solves a single MoMA problem", {
    # four grid requests

    zero_cnt <- 0

    for (av in alpha_v) {
        for (au in alpha_u) {
            for (lv in lambda_v) {
                for (lu in lambda_u) {
                    public_arglist_w_penalty <- c(
                        public_arglist_wo_penalty,
                        list(
                            alpha_u = au,
                            alpha_v = av,
                            lambda_u = lu,
                            lambda_v = lv
                        )
                    )

                    result <- do.call(
                        cpp_multirank_BIC_grid_search,
                        public_arglist_w_penalty
                    )

                    result_directcall <- do.call(
                        cpp_moma_multi_rank,
                        public_arglist_w_penalty
                    )

                    if (sum(result_directcall$u) == 0) {
                        zero_cnt <- zero_cnt + 1
                    }

                    expect_equal(result_directcall$u, get_5Dlist_elem(result, 1, 1, 1, 1)[[1]]$u$vector)
                    expect_equal(result_directcall$v, get_5Dlist_elem(result, 1, 1, 1, 1)[[1]]$v$vector)
                }
            }
        }
    }


    # Make sure penalty levels spread evenly from
    # those that zeros everything to those that are
    # all zeros
    expect_lte(zero_cnt / (n_lambda_u * n_lambda_v * n_alpha_u * n_alpha_v), 0.1)
})

test_that("`cpp_multirank_BIC_grid_search`: a naive case", {
    arglist_x_w_penalty <- modifyList(
        public_arglist_wo_penalty,
        list(
            X = matrix(1),
            Omega_u = matrix(0),
            Omega_v = matrix(0),
            alpha_u = 0,
            alpha_v = 0,
            lambda_u = 0,
            lambda_v = 0
        )
    )

    result <- do.call(
        cpp_multirank_BIC_grid_search,
        arglist_x_w_penalty
    )[[1]]

    expect_equal(result$u$vector, as.matrix(1))
    expect_equal(result$v$vector, as.matrix(1))
    expect_equal(result$d, 1)
})
