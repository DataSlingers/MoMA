context("Test R6 object")
set.seed(12)
test_that("SFPCA object: a naive case, 1x1 matrix", {
    a <- SFPCA$new(matrix(1), center = FALSE)
    result_to_be_tested <- a$grid_result[[1]]

    expect_equal(result_to_be_tested$X, matrix(1))

    expect_equal(result_to_be_tested$u$bic, -Inf)
    expect_equal(result_to_be_tested$u$lambda, 0)
    expect_equal(result_to_be_tested$u$alpha, 0)
    expect_equal(result_to_be_tested$u$vector, matrix(1))

    expect_equal(result_to_be_tested$v$bic, -Inf)
    expect_equal(result_to_be_tested$v$alpha, 0)
    expect_equal(result_to_be_tested$v$lambda, 0)
    expect_equal(result_to_be_tested$v$vector, matrix(1))

    expect_equal(result_to_be_tested$k, 0)
    expect_equal(result_to_be_tested$d, 1)
})

test_that("SFPCA object: correct arguments", {
    a <- SFPCA$new(matrix(runif(12), 3, 4), scale = FALSE, cen = FALSE)
    expect_equal(a$Omega_u, diag(3))
    expect_equal(a$Omega_v, diag(4))
    expect_equal(a$sc, NULL)

    a <- SFPCA$new(matrix(runif(12), 3, 4), alpha_u = 1)
    expect_equal(a$Omega_u, second_diff_mat(3))
    expect_equal(a$Omega_v, diag(4))

    a <- SFPCA$new(matrix(runif(12), 3, 4), u_sparsity = lasso())
    expect_equal(a$u_sparsity, lasso())
    expect_equal(a$v_sparsity, empty())

    expect_error(
        SFPCA$new(matrix(runif(12), 3, 4), selection_scheme_str = "bbba"),
        "Invalid selection_scheme_str bbba. It should be a four-char string containing only 'b' or 'g'. "
    )

    a <- SFPCA$new(matrix(runif(12), 3, 4), selection_scheme_str = "ggbb") # in order of alpha_u/v, lambda_u/v

    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 1)))

    expect_equal(dim(a$grid_result), c(1, 1, 1, 1, 1))
})



test_that("SFPCA object: as SVD", {
    X <- matrix(runif(12), 4, 3)
    expect_error(
        SFPCA$new(X, rank = "3.1", center = FALSE, scale = FALSE),
        "rank should be a legit positive integer"
    )
    expect_error(
        SFPCA$new(X, rank = 0, center = FALSE, scale = FALSE),
        "rank should be a legit positive integer"
    )
    expect_error(
        SFPCA$new(X, rank = -1.1, center = FALSE, scale = FALSE),
        "rank should be a legit positive integer"
    )
    expect_error(
        SFPCA$new(X, rank = -1, center = FALSE, scale = FALSE),
        "rank should be a legit positive integer"
    )
    expect_error(
        SFPCA$new(X, rank = 3.1, center = FALSE, scale = FALSE),
        "rank should be a legit positive integer"
    )

    # test no penalty case
    X <- matrix(runif(12), 4, 3)
    a <- SFPCA$new(X, rank = 3, center = FALSE, scale = FALSE)


    mysvd <- a$get_mat_by_index()

    svda <- svd(X, nu = 3, nv = 3)

    expect_equal(abs(mysvd$U), abs(svda$u), check.attributes = FALSE) # use `abs` to avoid sign difference
    expect_equal(abs(mysvd$V), abs(svda$v), check.attributes = FALSE)


    # test a matrix with dimnames
    rown <- paste0("row", seq(1:4))
    coln <- paste0("col", seq(1:3))
    dimnames(X) <- list(rown, coln)

    a <- SFPCA$new(X, rank = 3, center = FALSE, scale = FALSE)
    mysvd <- a$get_mat_by_index()

    expect_equal(rownames(mysvd$U), rown)
    expect_equal(rownames(mysvd$V), coln)


    # test default arguments
    expect_equal(a$get_mat_by_index(), a$get_mat_by_index())

    # test selection with BIC seach and grid search
    a <- SFPCA$new(X,
        rank = 3, center = FALSE, scale = FALSE,
        alpha_u = seq(0, 2, 0.2), selection_scheme_str = "bggg"
    )
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))

    expect_equal(dim(a$grid_result), c(1, 1, 1, 1, 3))
})

test_that("SFPCA object: print fucntion", {
    X <- matrix(runif(12), 4, 3)
    a <- SFPCA$new(X,
        rank = 3,
        alpha_u = c(1, 2),
        lambda_u = c(2, 3),
        alpha_v = c(3, 4),
        lambda_v = c(6, 7),
        selection_scheme_str = "bgbg"
    )
    print_message <- capture.output(print(a))

    expect_message <- c(
        "An <SFPCA> object containing solutions to the following settings",
        "rank = 3 ",
        "Penalty and selection:",
        "alpha_u: BIC search ",
        "[1] 1 2",
        "alpha_u: grid search ",
        "[1] 3 4",
        "lambda_u: BIC search ",
        "[1] 2 3",
        "lambda_v: grid search ",
        "[1] 6 7"
    )
    expect_equal(expect_message, print_message)
})

test_that("SFPCA object: left-project fucntion", {

    # unamed matrix
    X <- matrix(runif(17 * 8), 17, 8)
    a <- SFPCA$new(X,
        rank = 3,
        alpha_u = c(1),
        lambda_u = c(2),
        alpha_v = c(3),
        lambda_v = c(6),
        selection_scheme_str = "bbbb"
    )
    SFPCA$debug("left_project")
    expect_error(
        a$left_project(matrix(0, 4, 1)),
        "`newX` is incompatible with orignal data."
    )

    new_data <- matrix(runif(24), 3, 8)
    res <- a$left_project(new_data)

    V <- res$V
    # we verify the the projected data
    # satisfies the normal equation:
    # X^T X b = X^T y
    expect_equal(
        t(V) %*% V %*% t(res$proj_data),
        t(V) %*% t(res$scaled_data)
    )
})

test_that("SFPCA object: `fixed_list` functions as expected", {
    set.seed(113)
    X <- matrix(runif(17 * 8), 17, 8) * 10

    invaid_indices_error <-
        paste0(
            "Invalid index in SFPCA::get_mat_by_index. Do not specify indexes of parameters ",
            "i) that are chosen by BIC, or ",
            "ii) that are not specified during initialization of the SFCPA object, or ",
            "iii) that are scalars during initialization of the SFCPA object."
        )

    # case 1:
    # parameters that did not appear in initialization
    # should not appear in SFPCA::get_mat_by_index at all
    a <- moma_sfpca(X)
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 1),
        invaid_indices_error
    )
    expect_warning(
        a$get_mat_by_index(alphau = 1),
        paste0("extra argument ", sQuote("alphau"), " will be disregarded")
    )


    # case 2:
    # parameters that were scalars in initialization
    # should not appear in SFPCA::get_mat_by_index either
    a <- moma_sfpca(X, alpha_u = 1)
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 1),
        invaid_indices_error
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 2),
        invaid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 1),
        invaid_indices_error
    )


    # case 2:
    # parameters that were selected by BIC
    # should not appear in SFPCA::get_mat_by_index either
    a <- moma_sfpca(X,
        alpha_u = c(1, 2),
        alpha_v = c(1, 2, 3), # selected by BIC
        selection_scheme_str = "gbgg"
    )
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_no_error(
        a$get_mat_by_index(alpha_u = 2)
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 0),
        invaid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 0),
        invaid_indices_error
    )

    a <- moma_sfpca(X, alpha_u = c(1, 2))
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_no_error(
        a$get_mat_by_index(alpha_u = 2)
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 2, alpha_v = 0),
        invaid_indices_error
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 0),
        invaid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 0),
        invaid_indices_error
    )
})



test_that("SFPCA object wrappers: moma_spca", {
    X <- matrix(runif(17 * 8), 17, 8)

    # test inputs
    expect_error(
        moma_spca(X, lambda_v = c(1, 2), lambda_u = c(2, 3)),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X, lambda_v = c(1, 2), u_sparsity = empty()),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X, v_sparsity = empty(), u_sparsity = empty()),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X,
            v_sparsity = empty(), lambda_v = c(1, 2),
            lambda_u = c(0)
        ),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X, lambda_v = c(1, 2), lambda_u = c(0)),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_warning(moma_spca(X), "No sparsity is imposed!")


    # test selection schemes
    expect_error(
        moma_spca(X,
            lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
            selection_scheme_str = "gb"
        ),
        "`selection_scheme_str` should be either 'g' or 'b'"
    )

    expect_error(
        moma_spca(X,
            lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
            selection_scheme_str = "c"
        ),
        "`selection_scheme_str` should be either 'g' or 'b'"
    )


    expect_error(a <- moma_spca(X, lambda_u = c()),
        paste0(
            "All penalty levels (",
            sQuote("lambda_u"), ", ",
            sQuote("lambda_v"), ", ",
            sQuote("alpha_u"), ", ",
            sQuote("alpha_v"),
            ") must be numeric."
        ),
        fixed = TRUE
    )

    a <- moma_spca(X)
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X, selection_scheme_str = "b") # no effects
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X,
        lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
        selection_scheme_str = "g"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X,
        lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
        selection_scheme_str = "b"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))


    a <- moma_spca(X,
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "g"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X,
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "b"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))
})


test_that("SFPCA object wrappers: moma_twspca", {
    X <- matrix(runif(17 * 8), 17, 8)

    # test inputs
    expect_no_error(
        moma_twspca(X, lambda_v = c(1, 2), lambda_u = c(2, 3))
    )
    expect_no_error(
        moma_twspca(X, lambda_v = c(1, 2), u_sparsity = empty())
    )
    expect_no_error(
        moma_twspca(X, v_sparsity = empty(), u_sparsity = empty())
    )

    expect_no_error(
        moma_twspca(X,
            v_sparsity = empty(), lambda_v = c(1, 2),
            lambda_u = c(0)
        )
    )

    expect_warning(
        moma_twspca(X),
        "No sparsity is imposed!"
    )
    expect_warning(
        moma_twspca(X, lambda_u = seq(0, 2, 0.2), u_sparsity = lasso()),
        "Please use `moma_spca` if only one side is penalized"
    )


    # test selection schemes
    expect_error(
        moma_twspca(X,
            lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
            selection_scheme_str = "gbb"
        ),
        "`selection_scheme_str` should be of length two."
    )

    expect_error(
        moma_twspca(X,
            lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
            selection_scheme_str = "cc"
        ),
        "`selection_scheme_str` should consist of 'g' or 'b'"
    )

    expect_no_error(
        moma_twspca(X,
            lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
            selection_scheme_str = "bb"
        )
    )


    expect_warning(a <- moma_twspca(X))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_warning(a <- moma_twspca(X, selection_scheme_str = "bb"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 1)))
    expect_warning(a <- moma_twspca(X, selection_scheme_str = "gg"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    expect_warning(a <- moma_twspca(X, selection_scheme_str = "bg"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))
    expect_warning(a <- moma_twspca(X, selection_scheme_str = "gb"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))


    a <- moma_twspca(X,
        lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "bb"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 1)))

    a <- moma_twspca(X,
        lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "gg"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_twspca(X,
        lambda_u = seq(0, 2, 0.2), u_sparsity = lasso(),
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "gb"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))

    a <- moma_twspca(X,
        lambda_v = seq(0, 2, 0.2), v_sparsity = lasso(),
        selection_scheme_str = "bg"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))
})


test_that("SFPCA object wrappers: moma_fpca", {
    X <- matrix(runif(17 * 8), 17, 8)

    # test inputs
    expect_error(
        moma_fpca(X, alpha_v = c(1, 2), alpha_u = c(2, 3)),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_error(
        moma_fpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = lasso(),
            selection_scheme_str = "g"
        ),
        "Omega_u/v is not a matrix."
    )

    # test when no penalty matrix is provided
    a <- moma_fpca(X, alpha_v = c(1, 2))
    expect_equal(a$Omega_u, diag(17)) # Omega is set to identiy mat if u or v is unpenalzied
    expect_equal(a$Omega_v, second_diff_mat(8)) # the default penalty matrix

    expect_warning(
        a <- moma_fpca(X),
        "No smoothness is imposed!"
    )
    expect_equal(a$Omega_u, diag(17))
    expect_equal(a$Omega_v, diag(8))

    expect_error(
        moma_fpca(X,
            Omega_v = 2 * diag(8), alpha_v = c(1, 2),
            alpha_u = c(0)
        ),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_error(
        moma_fpca(X, alpha_v = c(1, 2), alpha_u = c(0)),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_warning(moma_fpca(X), "No smoothness is imposed!")


    # test selection schemes
    # error when nchar(selection_scheme_str) != 1
    expect_error(
        moma_fpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = second_diff_mat(17),
            selection_scheme_str = "gb"
        ),
        "`selection_scheme_str` should be either 'g' or 'b'"
    )

    expect_error(
        moma_fpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = second_diff_mat(8),
            selection_scheme_str = "c"
        ),
        "`selection_scheme_str` should be either 'g' or 'b'"
    )


    a <- moma_fpca(X)
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_fpca(X, selection_scheme_str = "b") # no effects
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_fpca(X,
        alpha_u = seq(0, 2, 0.2),
        selection_scheme_str = "g"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_fpca(X,
        alpha_u = seq(0, 2, 0.2),
        selection_scheme_str = "b"
    )
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))
    expect_equal(a$Omega_u, second_diff_mat(17))

    a <- moma_fpca(X,
        alpha_v = seq(0, 2, 0.2),
        selection_scheme_str = "g"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    expect_equal(a$Omega_v, second_diff_mat(8))

    a <- moma_fpca(X,
        alpha_v = seq(0, 2, 0.2),
        selection_scheme_str = "b"
    )
    expect_true(all(a$selection_scheme_list == c(0, 1, 0, 0)))
})


test_that("SFPCA object wrappers: moma_twfpca", {
    X <- matrix(runif(17 * 8), 17, 8)

    # test inputs
    expect_no_error(
        moma_twfpca(X, alpha_v = c(1, 2), alpha_u = c(2, 3))
    )
    expect_no_error(
        moma_twfpca(X, alpha_v = c(1, 2))
    )
    expect_no_error(
        moma_twfpca(X, Omega_v = diag(8), Omega_u = 2.1 * second_diff_mat(17))
    )

    # incompatible Omega
    expect_error(
        moma_twfpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = 2.1 * second_diff_mat(11),
            selection_scheme_str = "bb"
        ),
        "Omega shoud be a compatible matrix. It should be of 17x17, but is actually 11x11"
    )

    expect_no_error(
        moma_twfpca(X,
            Omega_v = diag(8), alpha_v = c(1, 2),
            alpha_u = c(0)
        )
    )

    expect_warning(
        moma_twfpca(X),
        "No smoothness is imposed!"
    )
    expect_warning(
        moma_twfpca(X, alpha_u = seq(0, 2, 0.2), Omega_u = 2.3 * second_diff_mat(17)),
        "Please use `moma_fpca` if only one side is penalized"
    )


    # test selection schemes
    expect_error(
        moma_twfpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = lasso(),
            selection_scheme_str = "gbb"
        ),
        "`selection_scheme_str` should be of length two."
    )

    expect_error(
        moma_twfpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = lasso(),
            selection_scheme_str = "cc"
        ),
        "`selection_scheme_str` should consist of 'g' or 'b'"
    )

    expect_no_error(
        moma_twfpca(X,
            alpha_u = seq(0, 2, 0.2), Omega_u = 2.1 * second_diff_mat(17),
            selection_scheme_str = "bb"
        )
    )



    expect_warning(a <- moma_twfpca(X))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_warning(a <- moma_twfpca(X, selection_scheme_str = "bb"))
    expect_true(all(a$selection_scheme_list == c(1, 1, 0, 0)))
    expect_warning(a <- moma_twfpca(X, selection_scheme_str = "gg"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    expect_warning(a <- moma_twfpca(X, selection_scheme_str = "bg"))
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))
    expect_warning(a <- moma_twfpca(X, selection_scheme_str = "gb"))
    expect_true(all(a$selection_scheme_list == c(0, 1, 0, 0)))


    a <- moma_twfpca(X,
        alpha_u = seq(0, 2, 0.2), Omega_u = 1.1 * second_diff_mat(17),
        alpha_v = seq(0, 2, 0.2), Omega_v = 1.2 * second_diff_mat(8),
        selection_scheme_str = "bb"
    )
    expect_true(all(a$selection_scheme_list == c(1, 1, 0, 0)))
    expect_equal(a$Omega_u, 1.1 * second_diff_mat(17))
    expect_equal(a$Omega_v, 1.2 * second_diff_mat(8))

    a <- moma_twfpca(X,
        alpha_u = seq(0, 2, 0.2), Omega_u = 1.1 * second_diff_mat(17),
        alpha_v = seq(0, 2, 0.2), Omega_v = 1.2 * second_diff_mat(8),
        selection_scheme_str = "gg"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    expect_equal(a$Omega_u, 1.1 * second_diff_mat(17))
    expect_equal(a$Omega_v, 1.2 * second_diff_mat(8))

    a <- moma_twfpca(X,
        alpha_u = seq(0, 2, 0.2), Omega_u = 1.1 * second_diff_mat(17),
        alpha_v = seq(0, 2, 0.2), Omega_v = 1.2 * second_diff_mat(8),
        selection_scheme_str = "gb"
    )
    expect_true(all(a$selection_scheme_list == c(0, 1, 0, 0)))

    a <- moma_twfpca(X,
        alpha_v = seq(0, 2, 0.2), Omega_v = 2.1 * second_diff_mat(8),
        selection_scheme_str = "bg"
    )
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))
})

test_that("SFPCA object: get_mat_by_index and left_project takes non-ingeters", {
    set.seed(113)
    X <- matrix(runif(17 * 8), 17, 8) * 10

    a <- moma_sfpca(X,
        alpha_u = c(1, 2), alpha_v = c(1.1, 2.1),
        lambda_u = c(1.3, 2.3), lambda_v = c(1.4, 2.4)
    )

    expect_no_error(
        a$get_mat_by_index(alpha_u = 2)
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 2.1),
        "Non-integer input in SFPCA::get_mat_by_index"
    )

    expect_no_error(
        a$left_project(X)
    )
    expect_error(
        a$left_project(X, alpha_u = 2.1),
        "Non-integer input in SFPCA::left_project"
    )
})

test_that("SFPCA object: interpolate", {
    set.seed(113)
    X <- matrix(runif(17 * 8), 17, 8) * 10

    a <- moma_sfpca(X,
        u_sparsity = lasso(), v_sparsity = lasso(),
        lambda_u = c(1, 2), lambda_v = c(1, 2),
        alpha_u = c(1, 2), alpha_v = c(1, 2),
        Omega_u = second_diff_mat(17), Omega_v = second_diff_mat(8),
        selection_scheme_str = "gggb"
    )
    expect_error(
        a$interpolate(),
        "R6 ojbect SFPCA do not support interpolation when BIC selection scheme has been used."
    )

    a <- moma_sfpca(X,
        u_sparsity = lasso(), v_sparsity = lasso(),
        lambda_u = 1.2, lambda_v = 1,
        alpha_u = 0.7, alpha_v = 0.3,
        Omega_u = second_diff_mat(17), Omega_v = second_diff_mat(8)
    )
    expect_error(
        a$interpolate(
            alpha_u = 1, exact = TRUE
        ),
        "SFPCA::interpolate does not support exact mode unless all parameters are specified."
    )

    SFPCA$debug("interpolate")
    internal_call <- a$interpolate(
        alpha_u = 1, alpha_v = 1,
        lambda_v = 1, lambda_u = 1,
        exact = TRUE
    )

    dir_call <- moma_svd(
        scale(X, center = a$center, scale = a$scale),
        u_sparsity = lasso(), v_sparsity = lasso(),
        alpha_u = 1, alpha_v = 1,
        lambda_v = 1, lambda_u = 1,
        Omega_u = second_diff_mat(17), Omega_v = second_diff_mat(8)
    )
    expect_equal(internal_call$U, dir_call$u)
    expect_equal(internal_call$V, dir_call$v)
})
