context("4D List")

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

result <- do.call(testnestedBIC,
                  c(arg_list,
                    list(bic_search_alpha_u=0,  #grid
                         bic_search_alpha_v=0,  #grid
                         bic_search_lambda_u=0, #grid
                         bic_search_lambda_v=0)))  #grid

test_that("Test 4D List attribute", {
    expect_true(inherits(result, "MoMA_4D_list"))
})


test_that("Passing wrong argument", {
    expect_error(get_4Dlist_elem(c(1),1,1,1,1),
                   paste0(sQuote("x"), " should be a ",sQuote("MoMA_4D_list")," object"))
})

test_that("Access all elements", {
    # dim(result) = 7 3 5 2
    # n_au=7; n_av=5; n_lu=3; n_lv=2; # mutually prime
    cnt = 1
    for(i in 1:n_au){
        for(j in 1:n_lu){
            for(k in 1:n_av){
                for(l in 1:n_lv){
                    expect_equal(get_4Dlist_elem(result,i,j,k,l)[[1]],
                                 result[[cnt]])
                    cnt = cnt + 1
                }
            }
        }
    }
})


test_that("Error when accessing broundary", {
    expect_no_error(get_4Dlist_elem(result,n_au,n_lu,n_av,n_lv))


    # NOTE: R index starts from 1
    expect_error(get_4Dlist_elem(result,n_au,n_lu,n_av,0),
                 "Invalid index \\(7,3,5,0\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,n_au,n_lu,0,n_lv),
                 "Invalid index \\(7,3,0,2\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,n_au,0,n_av,n_lv),
                 "Invalid index \\(7,0,5,2\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,0,n_lu,n_av,n_lv),
                 "Invalid index \\(0,3,5,2\\), dim = c\\(7, 3, 5, 2\\)")


    expect_error(get_4Dlist_elem(result,n_au,n_lu,n_av,n_lv+1),
                 "Invalid index \\(7,3,5,3\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,n_au,n_lu,n_av+1,n_lv),
                 "Invalid index \\(7,3,6,2\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,n_au,n_lu+1,n_av,n_lv),
                 "Invalid index \\(7,4,5,2\\), dim = c\\(7, 3, 5, 2\\)")

    expect_error(get_4Dlist_elem(result,n_au+1,n_lu,n_av,n_lv),
                 "Invalid index \\(8,3,5,2\\), dim = c\\(7, 3, 5, 2\\)")
})


