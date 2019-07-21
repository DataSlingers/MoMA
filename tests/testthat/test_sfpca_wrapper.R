context("Test R6 object")
test_that("SFPCA object: a naive case, 1x1 matrix", {
    set.seed(12)

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

test_that("SFPCA object: X contains strings", {
    set.seed(12)
    X <- matrix(seq(1:12), 4, 3)
    X[1, 1] <- "abc"
    expect_error(SFPCA$new(X), "X must contain numbers only")

    X <- matrix(seq(1:12), 4, 3)
    X[1, 1] <- Inf
    expect_error(SFPCA$new(X), "X must not have NaN, NA, or Inf")
})

test_that("SFPCA object: correct arguments", {
    set.seed(12)

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
    set.seed(12)

    X <- matrix(runif(12), 4, 3)
    expect_error(
        SFPCA$new(X, rank = "3.1", center = FALSE, scale = FALSE),
        "`rank` should be a positive integer smaller than the rank of the data matrix."
    )
    expect_error(
        SFPCA$new(X, rank = 0, center = FALSE, scale = FALSE),
        "`rank` should be a positive integer smaller than the rank of the data matrix."
    )
    expect_error(
        SFPCA$new(X, rank = -1.1, center = FALSE, scale = FALSE),
        "`rank` should be a positive integer smaller than the rank of the data matrix."
    )
    expect_error(
        SFPCA$new(X, rank = -1, center = FALSE, scale = FALSE),
        "`rank` should be a positive integer smaller than the rank of the data matrix."
    )
    expect_error(
        SFPCA$new(X, rank = 3.1, center = FALSE, scale = FALSE),
        "`rank` should be a positive integer smaller than the rank of the data matrix."
    )

    # test no penalty case
    X <- matrix(runif(12), 4, 3)
    a <- SFPCA$new(X, rank = 3, center = FALSE, scale = FALSE)


    mysvd <- a$get_mat_by_index()

    svda <- svd(X, nu = 3, nv = 3)

    expect_equal(abs(mysvd$U), abs(svda$u), check.attributes = FALSE) # use `abs` to avoid sign difference
    expect_equal(abs(mysvd$V), abs(svda$v), check.attributes = FALSE)
    expect_equal(mysvd$d, svda$d, check.attributes = FALSE)

    # test a matrix with dimnames
    rown <- paste0("row", seq(1:4))
    coln <- paste0("col", seq(1:3))
    dimnames(X) <- list(rown, coln)

    a <- SFPCA$new(X, rank = 3, center = FALSE, scale = FALSE)
    mysvd <- a$get_mat_by_index()

    expect_equal(rownames(mysvd$U), rown)
    expect_equal(rownames(mysvd$V), coln)


    # `a` is specified without
    # any parameters (all four parameters default
    # to 0 actually), so users must not specify any
    # parameters when calling `get_mat_by_index`.
    # A bit stringent though.
    expect_error(
        a$get_mat_by_index(alpha_u = 1),
        "Invalid index in SFPCA::get_mat_by_index."
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 1),
        "Invalid index in SFPCA::get_mat_by_index."
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 1),
        "Invalid index in SFPCA::get_mat_by_index."
    )
    expect_error(
        a$get_mat_by_index(lambda_v = 1),
        "Invalid index in SFPCA::get_mat_by_index."
    )
    expect_no_error(
        a$get_mat_by_index(),
        "Invalid index in SFPCA::get_mat_by_index."
    )

    # test selection with BIC seach and grid search
    a <- SFPCA$new(X,
        rank = 3, center = FALSE, scale = FALSE,
        alpha_u = seq(0, 2, 0.2), selection_scheme_str = "bggg"
    )
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))

    expect_equal(dim(a$grid_result), c(1, 1, 1, 1, 3))
})

test_that("SFPCA object: print fucntion", {
    set.seed(12)

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

    expected_message <- c(
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
    expect_equal(expected_message, print_message)
})

test_that("SFPCA object: left-project fucntion", {
    set.seed(12)

    X <- matrix(runif(17 * 8), 17, 8) * 10
    a <- SFPCA$new(X,
        rank = 3,
        alpha_u = c(1),
        lambda_u = c(2),
        alpha_v = c(3),
        lambda_v = c(6),
        selection_scheme_str = "bbbb"
    )
    expect_error(
        a$left_project(matrix(0, 4, 1)),
        "`newX` is incompatible with orignal data. It must have 8 columns."
    )

    new_data <- matrix(runif(24), 3, 8)
    res <- a$left_project(new_data)

    V <- res$V
    # We verify the the projected data
    # satisfies the normal equation:
    # X^T X b = X^T y.
    expect_equal(
        t(V) %*% V %*% t(res$proj_data),
        t(V) %*% t(res$scaled_data)
    )
})

test_that("SFPCA object: `fixed_list` functions as expected", {
    set.seed(113)
    X <- matrix(runif(17 * 8), 17, 8) * 10

    invalid_indices_error <-
        paste0(
            "Invalid index in SFPCA::get_mat_by_index. Do not specify indexes of parameters ",
            "i) that are chosen by BIC, or ",
            "ii) that are not specified during initialization of the SFCPA object, or ",
            "iii) that are scalars during initialization of the SFCPA object."
        )

    # case 1:
    # parameters that did not appear in initialization
    # should not appear in `SFPCA::get_mat_by_index` at all
    a <- moma_sfpca(X)
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 1),
        invalid_indices_error
    )
    # when an unused argument is given, or a typo.
    expect_warning(
        a$get_mat_by_index(alphau = 1),
        paste0("extra argument ", sQuote("alphau"), " will be disregarded")
    )


    # case 2:
    # parameters that were scalars in initialization
    # should not appear in SFPCA::get_mat_by_index either
    a <- moma_sfpca(X, u_smooth = moma_smoothness(alpha = 1))
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 1),
        invalid_indices_error
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 2),
        invalid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 1),
        invalid_indices_error
    )


    # case 3:
    # parameters that were selected by BIC
    # should not appear in SFPCA::get_mat_by_index either
    a <- moma_sfpca(X,
        u_smooth = moma_smoothness(alpha = c(1, 2)),
        v_smooth = moma_smoothness(alpha = c(1, 2), select_scheme = "b")
    )
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_no_error(
        a$get_mat_by_index(alpha_u = 2)
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 0),
        invalid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 0),
        invalid_indices_error
    )

    a <- moma_sfpca(X,
        u_smooth = moma_smoothness(alpha = c(1, 2))
    )
    expect_no_error(
        a$get_mat_by_index()
    )
    expect_no_error(
        a$get_mat_by_index(alpha_u = 2)
    )
    expect_error(
        a$get_mat_by_index(alpha_u = 2, alpha_v = 0),
        invalid_indices_error
    )
    expect_error(
        a$get_mat_by_index(alpha_v = 0),
        invalid_indices_error
    )
    expect_error(
        a$get_mat_by_index(lambda_u = 0),
        invalid_indices_error
    )
})

test_that("SFPCA object wrappers: moma_spca", {
    set.seed(12)

    X <- matrix(runif(17 * 8), 17, 8) * 10

    # test inputs
    expect_error(
        moma_spca(X,
            v_sparse = moma_empty(lambda = c(1, 2)),
            u_sparse = moma_empty(lambda = c(1, 2))
        ),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X,
            v_sparse = moma_empty(),
            u_sparse = moma_empty()
        ),
        "Please use `moma_twspca` if both sides are penalized"
    )
    expect_error(
        moma_spca(X,
            v_sparse = lasso()
        ),
        paste0(
            "Invalid argument: ",
            sQuote("u_sparse / v_sparse"),
            ". They should be of class `moma_sparsity_type`."
        )
    )


    expect_warning(moma_spca(X), "No sparsity is imposed!")


    # test selection schemes
    expect_error(
        moma_spca(X,
            u_smooth = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "c")
        ),
        "‘select_scheme’ is not valid: c. It should be either `g` or `b`.",
        fixed = TRUE
    )

    expect_error(
        moma_spca(X,
            u_smooth = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "gg")
        ),
        "‘select_scheme’ is not valid: gg. It should be either `g` or `b`.",
        fixed = TRUE
    )


    expect_error(
        moma_spca(X,
            u_smooth = moma_lasso(lambda = c())
        ),
        "‘lambda’ is not valid: NULL",
        fixed = TRUE
    )

    expect_warning(
        a <- moma_spca(X),
        "No sparsity is imposed!"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_no_warning(
        a <- moma_spca(X, u_sparse = moma_empty()),
        "No sparsity is imposed!"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X, u_sparse = moma_lasso(select_scheme = "g"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X, u_sparse = moma_lasso(select_scheme = "b"))
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))

    a <- moma_spca(X,
        v_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "g")
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_spca(X,
        v_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "b")
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))
})

test_that("SFPCA object wrappers: moma_twspca", {
    set.seed(12)

    X <- matrix(runif(17 * 8), 17, 8) * 10

    # test inputs
    expect_no_error(
        moma_twspca(X,
            v_sparse = moma_empty(lambda = c(1, 2)),
            u_sparse = moma_empty(lambda = c(2, 3))
        )
    )
    expect_no_error(
        moma_twspca(X,
            v_sparse = moma_empty(),
            u_sparse = moma_empty(lambda = c(2, 3))
        )
    )
    expect_no_error(
        moma_twspca(X,
            v_sparse = moma_empty(),
            u_sparse = moma_empty()
        )
    )

    expect_no_error(
        moma_twspca(X,
            v_sparse = moma_empty(lambda = c(0)),
            u_sparse = moma_empty()
        )
    )

    expect_warning(
        moma_twspca(X),
        "No sparsity is imposed!"
    )

    expect_warning(
        moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 2, 0.2))
        ),
        "Please use `moma_spca` if only one side is penalized"
    )


    # test selection schemes
    expect_warning(
        expect_error(
            moma_twspca(X,
                u_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "bb")
            ),
            "‘select_scheme’ is not valid: bb. It should be either `g` or `b`"
        ),
        "Please use `moma_spca` if only one side is penalized."
    )

    expect_warning(
        expect_error(
            moma_twspca(X,
                u_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "c")
            ),
            "‘select_scheme’ is not valid: c. It should be either `g` or `b`"
        ),
        "Please use `moma_spca` if only one side is penalized."
    )

    expect_warning(
        expect_no_error(
            moma_twspca(X,
                u_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "b")
            )
        ),
        "Please use `moma_spca` if only one side is penalized."
    )



    expect_warning(a <- moma_twspca(X))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "b"),
            v_sparse = moma_lasso(lambda = seq(0, 2, 0.2), select_scheme = "b")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 1)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2)),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2))
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))


    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b"),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2))
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2)),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))


    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b"),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 1)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "g"),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "g")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "g"),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 1)))

    expect_no_warning(
        a <- moma_twspca(X,
            u_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "b"),
            v_sparse = moma_lasso(lambda = seq(0, 1, 0.2), select_scheme = "g")
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 1, 0)))
})

test_that("SFPCA object wrappers: moma_fpca", {
    set.seed(12)

    X <- matrix(runif(17 * 8), 17, 8) * 10

    # test inputs
    expect_error(
        moma_fpca(X,
            u_smooth = moma_smoothness(alpha = c(1, 2)),
            v_smooth = moma_smoothness(alpha = c(1, 2))
        ),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_error(
        moma_fpca(X,
            u_smooth = moma_smoothness(
                Omega = lasso(),
                alpha = c(1, 2)
            )
        ),
        "Omega_u/v is not a matrix."
    )

    # test when no penalty matrix is provided
    a <- moma_fpca(X,
        v_smooth = moma_smoothness(alpha = c(1, 2))
    )
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
            v_smooth = moma_smoothness(Omega = 2 * diag(8), alpha = c(1, 2)),
            u_smooth = moma_smoothness(alpha = c(0))
        ),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_error(
        moma_fpca(X,
            v_smooth = moma_smoothness(alpha = c(1, 2)),
            u_smooth = moma_smoothness(alpha = c(0))
        ),
        "Please use `moma_twfpca` if both sides are penalized"
    )
    expect_warning(moma_fpca(X), "No smoothness is imposed!")


    # test selection schemes
    # error when nchar(selection_scheme_str) != 1
    expect_error(
        moma_fpca(X,
            u_smooth = moma_smoothness(
                Omega = second_diff_mat(17),
                alpha = seq(0, 2, 0.2),
                select_scheme = "bg"
            )
        ),
        "‘select_scheme’ is not valid: bg. It should be either `g` or `b`",
        fixed = TRUE
    )

    expect_error(
        moma_fpca(X,
            u_smooth = moma_smoothness(
                Omega = second_diff_mat(17),
                alpha = seq(0, 2, 0.2),
                select_scheme = "c"
            )
        ),
        "‘select_scheme’ is not valid: c. It should be either `g` or `b`",
        fixed = TRUE
    )


    expect_warning(
        a <- moma_fpca(X),
        "No smoothness is imposed!"
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))


    a <- moma_fpca(X, u_smooth = moma_smoothness(alpha = seq(0, 2, 0.2)))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    a <- moma_fpca(X, u_smooth = moma_smoothness(
        alpha = seq(0, 2, 0.2),
        select_scheme = "b"
    ))

    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))
    expect_equal(a$Omega_u, second_diff_mat(17))

    a <- moma_fpca(X,
        v_smooth = moma_smoothness(alpha = seq(0, 2, 0.2))
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    expect_equal(a$Omega_v, second_diff_mat(8))
    expect_equal(a$Omega_u, diag(17))


    a <- moma_fpca(X,
        v_smooth = moma_smoothness(
            alpha = seq(0, 2, 0.2),
            select_scheme = "b"
        )
    )
    expect_true(all(a$selection_scheme_list == c(0, 1, 0, 0)))
})

test_that("SFPCA object wrappers: moma_twfpca", {
    set.seed(12)

    X <- matrix(runif(17 * 8), 17, 8) * 10

    # test inputs
    expect_no_error(
        moma_twfpca(X,
            v_smooth = moma_smoothness(alpha = c(1, 2)),
            u_smooth = moma_smoothness(alpha = c(2, 3))
        )
    )

    expect_warning(
        expect_no_error(
            moma_twfpca(X,
                u_smooth = moma_smoothness(alpha = c(2, 3))
            )
        ),
        "Please use `moma_fpca` if only one side is penalized"
    )

    # incompatible Omega
    expect_warning(
        expect_error(
            moma_twfpca(X,
                u_smooth = moma_smoothness(Omega = 2.1 * second_diff_mat(11), alpha = seq(0, 2, 0.2))
            ),
            "Omega shoud be a compatible matrix. It should be of 17x17, but is actually 11x11"
        ),
        "Please use `moma_fpca` if only one side is penalized"
    )


    expect_no_error(
        moma_twfpca(X,
            v_smooth = moma_smoothness(Omega = diag(8), alpha = c(1, 2)),
            u_smooth = moma_smoothness(alpha = c(0))
        )
    )

    expect_warning(
        moma_twfpca(X),
        "No smoothness is imposed!"
    )
    expect_warning(
        moma_twfpca(X,
            u_smooth = moma_smoothness(Omega = 2.3 * second_diff_mat(17)),
            alpha = seq(0, 2, 0.2)
        ),
        "Please use `moma_fpca` if only one side is penalized"
    )


    expect_warning(a <- moma_twfpca(X))
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))

    expect_no_error(
        a <- moma_twfpca(X,
            v_smooth = moma_smoothness(
                Omega = diag(8),
                alpha = c(1, 2)
            ),
            u_smooth = moma_smoothness(
                Omega = second_diff_mat(8),
                alpha = c(0)
            ) # if alpha = 0, Omega will be overwritten,
        )
    )

    expect_error(
        a <- moma_twfpca(X,
            v_smooth = moma_smoothness(Omega = diag(8), alpha = c(1, 2)),
            u_smooth = moma_smoothness(Omega = second_diff_mat(8), alpha = c(1))
        ),
        "Omega shoud be a compatible matrix. It should be of 17x17, but is actually 8x8 "
    )

    a <- moma_twfpca(X,
        u_smooth = moma_smoothness(select_scheme = "g"),
        v_smooth = moma_smoothness(select_scheme = "b")
    )
    expect_true(all(a$selection_scheme_list == c(0, 1, 0, 0)))
    a <- moma_twfpca(X,
        u_smooth = moma_smoothness(select_scheme = "b"),
        v_smooth = moma_smoothness(select_scheme = "b")
    )
    expect_true(all(a$selection_scheme_list == c(1, 1, 0, 0)))
    a <- moma_twfpca(X,
        u_smooth = moma_smoothness(select_scheme = "g"),
        v_smooth = moma_smoothness(select_scheme = "g")
    )
    expect_true(all(a$selection_scheme_list == c(0, 0, 0, 0)))
    a <- moma_twfpca(X,
        u_smooth = moma_smoothness(select_scheme = "b"),
        v_smooth = moma_smoothness(select_scheme = "g")
    )
    expect_true(all(a$selection_scheme_list == c(1, 0, 0, 0)))
})

# test_that("SFPCA object: get_mat_by_index and left_project takes non-ingeters", {
#     set.seed(113)
#     X <- matrix(runif(17 * 8), 17, 8) * 10

#     a <- moma_sfpca(X,
#         alpha_u = c(1, 2), alpha_v = c(1.1, 2.1),
#         lambda_u = c(1.3, 2.3), lambda_v = c(1.4, 2.4)
#     )

#     expect_no_error(
#         a$get_mat_by_index(alpha_u = 2)
#     )
#     expect_error(
#         a$get_mat_by_index(alpha_u = 2.1),
#         "Non-integer input in SFPCA::get_mat_by_index"
#     )

#     expect_no_error(
#         a$left_project(X)
#     )
#     expect_error(
#         a$left_project(X, alpha_u = 2.1),
#         "Non-integer input in SFPCA::left_project"
#     )
# })

# test_that("SFPCA object: interpolate, exact mode", {
#     set.seed(113)
#     X <- matrix(runif(17 * 8), 17, 8) * 10

#     # Once BIC is used, `SFPCA::interpolate`
#     # cannot be called at all.
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         lambda_u = c(1, 2), lambda_v = c(1, 2),
#         alpha_u = c(1, 2),
#         alpha_v = c(1, 2), # selected by BIC
#         Omega_u = second_diff_mat(17), Omega_v = second_diff_mat(8),
#         selection_scheme_str = "gggb"
#     )
#     expect_error(
#         a$interpolate(),
#         "R6 object SFPCA do not support interpolation when BIC selection scheme has been used."
#     )
#     expect_error(
#         a$interpolate(
#             lambda_u = 1.1, lambda_v = 1.3,
#             alpha_u = 1.4
#         ),
#         "R6 object SFPCA do not support interpolation when BIC selection scheme has been used."
#     )

#     # case 1: `SFPCA::interpolate` cannot be used at all if all
#     # parameters are specified as scalars
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         lambda_u = 1.2, lambda_v = 1,
#         alpha_u = 0.7, alpha_v = 0.3,
#         Omega_u = second_diff_mat(17), Omega_v = second_diff_mat(8)
#     )
#     # this returns `a$get_mat_by_id()`
#     expect_no_error(
#         a$interpolate(
#             exact = TRUE
#         )
#     )
#     expect_error(
#         a$interpolate(
#             exact = FALSE
#         ),
#         "SFPCA::interpolate only supports one-sided interpolation"
#     )
#     expect_error(
#         a$interpolate(
#             alpha_u = 1, exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: alpha_u."
#     )
#     expect_error(
#         a$interpolate(
#             alpha_v = 1, exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: alpha_v."
#     )
#     expect_error(
#         a$interpolate(
#             lambda_u = 1, exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: lambda_u."
#     )
#     expect_error(
#         a$interpolate(
#             lambda_v = 1, exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: lambda_v."
#     )


#     # case 2: interpolation on v side, both alpha and lambda are vectors
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         alpha_v = seq(0.1, 1, 0.15), alpha_u = 0.32,
#         Omega_v = second_diff_mat(8), Omega_u = second_diff_mat(17),
#         lambda_v = seq(0.1, 1, 0.2), lambda_u = 2.1
#     )
#     # incomplete arguments
#     expect_error(
#         a$interpolate(
#             lambda_v = 1, exact = TRUE
#         ),
#         "Please spesify the following argument(s): alpha_v.",
#         fixed = TRUE
#     )
#     # extra arguments
#     expect_error(
#         a$interpolate(
#             lambda_v = 1, alpha_v = 1, alpha_u = 1,
#             exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: alpha_u."
#     )
#     # alpha_v too large
#     expect_error(
#         a$interpolate(
#             lambda_v = 0.21, alpha_v = 1
#         ),
#         "Invalid range: alpha_v."
#     )

#     # case 3: interpolation on v side, only lambda is a vector
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         alpha_v = 0.12, alpha_u = 0.32,
#         Omega_v = second_diff_mat(8), Omega_u = second_diff_mat(17),
#         lambda_v = seq(0.1, 1, 0.2), lambda_u = 2.1
#     )
#     # correct arguments
#     expect_no_error(
#         a$interpolate(
#             lambda_v = 1, exact = TRUE
#         )
#     )
#     # extra arguments
#     expect_error(
#         a$interpolate(
#             lambda_v = 1, alpha_v = 1, alpha_u = 1,
#             exact = TRUE
#         ),
#         "Invalid index in SFPCA::interpolate: alpha_u, alpha_v."
#     )
#     # alpha_v must not be specified
#     expect_error(
#         a$interpolate(
#             lambda_v = 0.21, alpha_v = 0.09
#         ),
#         "Invalid index in SFPCA::interpolate: alpha_v."
#     )
# })

# test_that("SFPCA object: interpolate, inexact mode", {
#     set.seed(113)
#     X <- matrix(runif(17 * 8), 17, 8) * 10

#     # interpolation on v side
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         alpha_v = seq(0.1, 1, 0.15), alpha_u = 0.32,
#         Omega_v = second_diff_mat(8), Omega_u = second_diff_mat(17),
#         lambda_v = seq(0.1, 1, 0.2), lambda_u = 2.1
#     )
#     expect_no_error(
#         a$interpolate(alpha_v = 0.23, lambda_v = 0.121)
#     )

#     # error because both alpha_v and lambda_v should be specified
#     expect_error(
#         a$interpolate(alpha_v = 0.23),
#         "Please spesify the following argument(s): lambda_v.",
#         fixed = TRUE
#     )
#     expect_error(
#         a$interpolate(),
#         "Please spesify the following argument(s): alpha_v, lambda_v.",
#         fixed = TRUE
#     )

#     # error because alpha_u is a scalar during initialization
#     expect_error(
#         a$interpolate(alpha_u = 0.2323),
#         "Invalid index in SFPCA::interpolate: alpha_u"
#     )

#     # error because alpha_u should not be specified
#     expect_error(
#         a$interpolate(alpha_v = 0.23, lambda_v = 0.121, alpha_u = 1),
#         "Invalid index in SFPCA::interpolate: alpha_u"
#     )
#     expect_error(
#         a$interpolate(alpha_v = 0.23, lambda_v = 0.121, alpha_u = 1, lambda_u = 1.3),
#         "Invalid index in SFPCA::interpolate: alpha_u, lambda_u."
#     )

#     # unsorted alpha_v
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         alpha_v = c(seq(0.1, 1, 0.15), 0.02),
#         alpha_u = 0.32,
#         Omega_v = second_diff_mat(8), Omega_u = second_diff_mat(17),
#         lambda_v = seq(0.1, 1, 0.2), lambda_u = 2.1
#     )

#     expect_error(
#         a$interpolate(alpha_v = 0.23, lambda_v = 0.121),
#         "Penalty levels not sorted"
#     )

#     # interpolate on both sides is not supported
#     a <- moma_sfpca(X,
#         u_sparsity = lasso(), v_sparsity = lasso(),
#         alpha_v = seq(0.1, 1, 0.15), alpha_u = seq(0.1, 1, 0.15),
#         Omega_v = second_diff_mat(8), Omega_u = second_diff_mat(17),
#         lambda_v = seq(0.1, 1.4, 0.4), lambda_u = seq(0.1, 1.3, 0.45)
#     )
#     # trying to do two-sided interpolation
#     # but not supported now.
#     expect_error(
#         a$interpolate(
#             alpha_v = 0.23, lambda_v = 0.121,
#             alpha_u = 0.1, lambda_u = 0.3
#         ),
#         "SFPCA::interpolate only supports one-sided interpolation"
#     )
# })
