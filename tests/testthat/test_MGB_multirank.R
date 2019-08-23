context("`cpp_multirank_BIC_grid_search` as multirank")

set.seed(123)
X <- matrix(runif(12), 3, 4) * 10

# They are mutually prime. This is useful for testing
# the size of the 5D list
n_alpha_u <- 2
n_alpha_v <- 3
n_lambda_u <- 5
n_lambda_v <- 7

# Penalties should not be too high. Otherwise
# u and v becomes vectors of zeros, and defaltion
# is thus trivial.
alpha_u <- seq(0, 1, length.out = n_alpha_u)
alpha_v <- seq(0, 1, length.out = n_alpha_v)
lambda_u <- seq(0, 1, length.out = n_lambda_u)
lambda_v <- seq(0, 1, length.out = n_lambda_v)

# public_arglist_wo_rank_and_selection does not specify two things: rank
# and selection strategy for each parameter. They
# are specified in each test case.
public_arglist_wo_rank_and_selection <- c(
    list(
        X = X,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(3), Omega_v = second_diff_mat(4),
        lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = add_default_prox_args(lasso()), prox_arg_list_v = add_default_prox_args(lasso())
    ),
    moma_pg_settings()
)

# Tests for multirank
test_that("cpp_multirank_BIC_grid_search receives rank <= 0", {
    # case 1: rank=0
    arglist_w_rank <- c(
        public_arglist_wo_rank_and_selection,
        list(
            rank = 0
        )
    )

    expect_error(
        do.call(cpp_multirank_BIC_grid_search, arglist_w_rank),
        "rank in MoMA::grid_BIC_mix should >= 1"
    )

    # case 2: rank=-1
    arglist_w_rank$rank <- -1
    expect_error(
        do.call(cpp_multirank_BIC_grid_search, arglist_w_rank),
        "rank in MoMA::grid_BIC_mix should >= 1"
    )
})

test_that("cpp_multirank_BIC_grid_search returns multi-rank solution", {
    # four grid requests
    rank <- 3

    public_arglist_w_rank <-
        modifyList(
            public_arglist_wo_rank_and_selection,
            list(rank = rank)
        )

    result <- do.call(
        cpp_multirank_BIC_grid_search,
        public_arglist_w_rank
    )

    for (i in 1:n_alpha_u) {
        for (j in 1:n_alpha_v) {
            for (n in 1:n_lambda_u) {
                for (m in 1:n_lambda_v) {
                    arglist_unary_penalty_wo_selection <-
                        modifyList(
                            public_arglist_wo_rank_and_selection,
                            list(
                                alpha_u = alpha_u[i],
                                alpha_v = alpha_v[j],
                                lambda_u = lambda_u[n],
                                lambda_v = lambda_v[m],
                                rank = rank
                            )
                        )

                    # call cpp_moma_multi_rank
                    result_directcall <- do.call(
                        cpp_moma_multi_rank,
                        arglist_unary_penalty_wo_selection
                    )

                    for (rank_i in 1:rank) {
                        res <- get_5Dlist_elem(result, alpha_u_i = i, lambda_u_i = n, alpha_v_i = j, lambda_v_i = m, rank_i)[[1]]
                        expect_equal(abs(matrix(result_directcall$u[, rank_i])), abs(res$u$vector))
                    }
                }
            }
        }
    }
})
