context("5D List")

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
alpha_u <- seq(0, 0.1, length.out = n_alpha_u)
alpha_v <- seq(0, 0.2, length.out = n_alpha_v)
lambda_u <- seq(0, 0.3, length.out = n_lambda_u)
lambda_v <- seq(0, 0.4, length.out = n_lambda_v)


public_arg_list <- c(
    list(
        X = X,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = second_diff_mat(3), Omega_v = second_diff_mat(4),
        lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = add_default_prox_args(lasso()), prox_arg_list_v = add_default_prox_args(empty()),
        rank = rank
    ),
    moma_pg_settings()
)

# Generate a 5-D list by calling `cpp_multirank_BIC_grid_search`
five_D_list_instance <- do.call(
    cpp_multirank_BIC_grid_search,
    c(
        public_arg_list,
        list(
            select_scheme_alpha_u = 0, # grid
            select_scheme_alpha_v = 0, # grid
            select_scheme_lambda_u = 0, # grid
            select_scheme_lambda_v = 0 # grid
        )
    )
)

test_that("Test 5D List attribute", {
    expect_true(inherits(five_D_list_instance, "MoMA_5D_list"))
})


test_that("Error on receiving non-MoMA_5D_list object", {
    expect_error(
        get_5Dlist_elem(c(1), 1, 1, 1, 1),
        paste0(sQuote("x"), " should be a ", sQuote("MoMA_5D_list"), " object")
    )
})


test_that("Access all elements", {
    # dim(five_D_list_instance) = 7 3 5 2
    # n_alpha_u=7; n_alpha_v=5; n_lambda_u=3; n_lambda_v=2; # mutually prime
    cnt <- 1
    for (i in 1:n_alpha_u) {
        for (j in 1:n_lambda_u) {
            for (k in 1:n_alpha_v) {
                for (l in 1:n_lambda_v) {
                    for (rk in 1:rank) {
                        expect_true(!(is.null(five_D_list_instance[[cnt]])))
                        expect_equal(
                            get_5Dlist_elem(five_D_list_instance, i, j, k, l, rk)[[1]],
                            five_D_list_instance[[cnt]]
                        )
                        cnt <- cnt + 1
                    }
                }
            }
        }
    }
})

test_that("No error when accessing the broundary", {
    expect_no_error(get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v, rank))
})

test_that("Error when just crossing the broundary", {
    # NOTE: R index starts from 1

    # Lower boundary
    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u, n_alpha_v, 0),
        "Invalid index \\(7,3,5,0,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u, 0, n_lambda_v),
        "Invalid index \\(7,3,0,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, 0, n_alpha_v, n_lambda_v),
        "Invalid index \\(7,0,5,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, 0, n_lambda_u, n_alpha_v, n_lambda_v),
        "Invalid index \\(0,3,5,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    # Upper boundary
    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u, n_alpha_v, n_lambda_v + 1),
        "Invalid index \\(7,3,5,3,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u, n_alpha_v + 1, n_lambda_v),
        "Invalid index \\(7,3,6,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u, n_lambda_u + 1, n_alpha_v, n_lambda_v),
        "Invalid index \\(7,4,5,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )

    expect_error(
        get_5Dlist_elem(five_D_list_instance, n_alpha_u + 1, n_lambda_u, n_alpha_v, n_lambda_v),
        "Invalid index \\(8,3,5,2,1\\), dim = c\\(7, 3, 5, 2, 3\\)"
    )
})
