context("UX / Argument parsing tests")
test_that("Test for arguments names", {
    # This test makes sure the lasso(), mcp(), etc.,
    # return arguments with
    # correct names.

    correct_args <- names(formals(sfpca))

    # Collect all the arguments
    test_args <- list()

    for (fun in c(
        lasso,
        scad,
        mcp,
        fusedlasso,
        l1tf,
        slope
    )) {
        test_args <- c(test_args, names(fun()))
    }

    args <- names(grplasso(g = factor(rep(1:10))))
    test_args <- c(test_args, args)

    args <- names(cluster(w = matrix(1)))
    test_args <- c(test_args, args)

    # Test prox arguments
    for (arg in test_args) {
        expect_true(paste0(arg, "_u") %in% correct_args)
        expect_true(paste0(arg, "_v") %in% correct_args)
    }

    # Test PG loop arguments
    for (arg in names(moma_pg_settings())) {
        expect_true(arg %in% correct_args)
    }
})

test_that("Prompt errors for wrong prox arguments", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))


    # Wrong non-convexity arguments
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = scad(gamma = 1),
            lambda_u = 3
        ),
        paste0(
            "Non-convexity parameter of SCAD (",
            sQuote("gamma"),
            ") must be larger than 2."
        ),
        fixed = TRUE
    )
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = mcp(gamma = 0.9),
            lambda_u = 3
        ),
        paste0(
            "Non-convexity parameter of MCP (",
            sQuote("gamma"),
            ") must be larger than 1."
        ),
        fixed = TRUE
    )


    # Wrong grouping dimension in group lasso
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = grplasso(g = factor(1)),
            lambda_u = 3
        ),
        "Wrong dimension: length(group) != dim(x).",
        fixed = TRUE
    )
    expect_error(
        grplasso(g = matrix(1)),
        "Please provide a vector as an indicator of grouping. (Called from grplasso)",
        fixed = TRUE
    )


    # Wrong weight matrix dimension in cluster penalty
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = cluster(w = matrix(1)),
            lambda_u = 3
        ),
        "Wrong dimension: dim(weight matrix) != dim(x).",
        fixed = TRUE
    )


    # Omega has wrong dimension
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            Omega_u = matrix(c(1, 2), 1, 2),
            alpha_u = 2
        ),
        "Omega shoud be a square matrix: nrows = 1, ncols = 2 (Called from check_omega)",
        fixed = TRUE
    )
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            Omega_u = matrix(c(1), 1),
            alpha_u = 2
        ),
        "Omega shoud be a compatible matrix. It should be of 3x3, but is actually 1x1 (Called from check_omega)",
        fixed = TRUE
    )


    # Prompt errors when users require rank-k svd and cross validation
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            lambda_u = c(1, 2, 3),
            k = 2
        ),
        "We don't support a range of parameters in finding a rank-k svd (Called from moma_svd)",
        fixed = TRUE
    )
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            lambda_u = c(1, 2, 3),
            lambda_v = seq(10),
            alpha_u = seq(10)
        ),
        "We only allow changing two parameters.",
        fixed = TRUE
    )

    expect_no_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            lambda_u = c(1, 2, 3),
            lambda_v = seq(10), select = "grid"
        )
    )

    expect_no_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            lambda_u = c(1, 2, 3),
            lambda_v = seq(10), select = "nested"
        )
    )
})


test_that("Correct prox match", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))

    expect_output(
        moma_svd(matrix(runif(12), 3, 4)),
        "Initializing null proximal operator object"
    )

    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = lasso(),
            lambda_u = 3
        ),
        "Initializing Lasso proximal operator object"
    )

    # SLOPE
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = slope(),
            lambda_u = 3
        ),
        "P_u SLOPE"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = slope(),
            lambda_u = 3
        ),
        "Initializing SLOPE proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            v_sparsity = slope(),
            lambda_u = 3
        ),
        "P_v SLOPE"
    )
    expect_output(
        moma_svd(matrix(runif(12), 3, 4),
            v_sparsity = slope()
        ),
        "Initializing SLOPE proximal operator object"
    )

    # lasso
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = lasso(),
            lambda_u = 3
        ),
        "Initializing Lasso proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = lasso(non_negative = TRUE),
            lambda_u = 3
        ),
        "Initializing non-negative Lasso proximal operator object"
    )


    # scad
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = scad(),
            lambda_u = 3
        ),
        "Initializing SCAD proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = scad(non_negative = TRUE),
            lambda_u = 3
        ),
        "Initializing non-negative SCAD proximal operator object"
    )


    # mcp
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = mcp(),
            lambda_u = 3
        ),
        "Initializing MCP proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = mcp(non_negative = TRUE),
            lambda_u = 3
        ),
        "Initializing non-negative MCP proximal operator object"
    )


    # group
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = grplasso(g = factor(seq(3))),
            lambda_u = 3
        ),
        "Initializing group lasso proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = grplasso(g = factor(seq(3)), non_negative = TRUE),
            lambda_u = 3
        ),
        "Initializing non-negative group lasso proximal operator object"
    )


    # L1 linear trend filtering
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = l1tf(),
            lambda_u = 3
        ),
        "Initializing a L1 linear trend filtering proximal operator object of degree 1"
    )
    expect_output(
        moma_svd(
            matrix(runif(100), 10, 10),
            u_sparsity = l1tf(l1tf_k = 2),
            lambda_u = 3
        ),
        "Initializing a L1 linear trend filtering proximal operator object of degree 2"
    )
    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = l1tf(l1tf_k = 2),
            lambda_u = 3
        ),
        "A difference matrix should have more columns."
    )

    # sparse fused lasso
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = spfusedlasso(lambda2 = 3),
            lambda_u = 3
        ),
        "Initializing a sparse fused lasso proximal operator object"
    )
    # fused lasso
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(),
            lambda_u = 3
        ),
        "Initializing a ordered fusion lasso proximal operator object"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(),
            lambda_u = 3
        ),
        "P_u ORDEREDFUSED P_v NONE"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            v_sparsity = fusedlasso(),
            lambda_u = 3
        ),
        "P_u NONE P_v ORDEREDFUSED"
    )

    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(algo = "dp"),
            lambda_u = 3
        ),
        "P_u ORDEREDFUSEDDP P_v NONE"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(algo = "dp"),
            v_sparsity = fusedlasso(),
            lambda_u = 3
        ),
        "P_u ORDEREDFUSEDDP P_v ORDEREDFUSED"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(algo = "dp"),
            lambda_u = 3
        ),
        "P_u ORDEREDFUSEDDP P_v NONE"
    )

    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = fusedlasso(algo = "dp"),
            lambda_u = 3
        ),
        "Initializing a ordered fusion lasso proximal operator object \\(DP\\)"
    )
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            v_sparsity = fusedlasso(algo = "dp"),
            lambda_u = 3
        ),
        "Initializing a ordered fusion lasso proximal operator object \\(DP\\)"
    )

    # cluster penalty
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = cluster(w = diag(3)),
            lambda_u = 3
        ),
        "Initializing a fusion lasso proximal operator object"
    )

    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = cluster(w = diag(3), ADMM = TRUE),
            lambda_u = 3
        ),
        "Running ADMM"
    )

    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            u_sparsity = cluster(w = diag(3)),
            lambda_u = 3
        ),
        "Running AMA"
    )
})

test_that("Correct match for PG loop settings", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")


    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "ista")),
        "Initializing a ISTA solver"
    )
    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "fista")),
        "Initializing a FISTA solver"
    )
    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "onestepista")),
        "Initializing an one-step ISTA solver"
    )


    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "ista")),
        "Releasing a ISTA object"
    )
    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "fista")),
        "Releasing a FISTA object"
    )
    expect_output(
        moma_svd(matrix(runif(12), 3, 4), pg_settings = moma_pg_settings(solver = "onestepista")),
        "Releasing a OneStepISTA object"
    )

    # Test default PG loop setting
    expect_output(
        moma_svd(matrix(runif(12), 3, 4)),
        "EPS 1e-10 MAX_ITER 1000 EPS_inner 1e-10 MAX_ITER_inner 100000 solver ISTA"
    )

    # Test pg_settings() passes correct values to C++ side
    expect_output(
        moma_svd(
            matrix(runif(12), 3, 4),
            pg_settings = moma_pg_settings(
                EPS = 1.212312e-5,
                MAX_ITER = 1.2957e+7,
                EPS_inner = 1.987e-6,
                MAX_ITER_inner = 98728376
            )
        ),
        "EPS 1.21231e-05 MAX_ITER 12957000 EPS_inner 1.987e-06 MAX_ITER_inner 98728376"
    )

    expect_error(
        moma_svd(
            matrix(runif(12), 3, 4),
            pg_settings = c(
                EPS = 1.212312e-5,
                MAX_ITER = 1.2957e+7,
                EPS_inner = 1.987e-6,
                MAX_ITER_inner = 98728376
            )
        ),
        paste0("pg_settings penalty should be of class ", sQuote("moma_pg_settings"))
    )

    on.exit(MoMA::moma_logger_level(old_logger_level))
})

test_that("Data matrix must be complete", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))

    X <- matrix(runif(12), 3, 4)
    X[2, 1] <- NA
    expect_error(moma_svd(X = X),
        "X must not have NaN, NA, or Inf. (Called from moma_svd)",
        fixed = TRUE
    )

    X <- matrix(runif(12), 3, 4)
    X[3, 2] <- Inf
    expect_error(moma_svd(X = X),
        "X must not have NaN, NA, or Inf. (Called from moma_svd)",
        fixed = TRUE
    )

    X <- matrix(runif(12), 3, 4)
    X[1, 4] <- NaN
    expect_error(moma_svd(X = X),
        "X must not have NaN, NA, or Inf. (Called from moma_svd)",
        fixed = TRUE
    )
})

test_that("Negative penalty", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))


    # Negative penalty
    set.seed(112)
    X <- matrix(runif(12), 3, 4)
    expect_error(
        moma_svd(X = X, lambda_u = c(0, 1, 2, 3, 4, -1)),
        paste0(
            "All penalty levels (",
            sQuote("lambda_u"),
            ", ",
            sQuote("lambda_v"),
            ", ",
            sQuote("alpha_u"),
            ", ",
            sQuote("alpha_v"),
            ") must be non-negative numeric. "
        ),
        fixed = TRUE
    )


    expect_error(
        moma_svd(X = X, lambda_v = c(0, 1, 2, 3, 4, -1)),
        paste0(
            "All penalty levels (",
            sQuote("lambda_u"),
            ", ",
            sQuote("lambda_v"),
            ", ",
            sQuote("alpha_u"),
            ", ",
            sQuote("alpha_v"),
            ") must be non-negative numeric. "
        ),
        fixed = TRUE
    )

    # Prompt error when penalty contains Infty
    expect_error(
        moma_svd(X = X, lambda_v = c(1:3, Inf)),
        paste0(
            "All penalty levels (",
            sQuote("lambda_u"),
            ", ",
            sQuote("lambda_v"),
            ", ",
            sQuote("alpha_u"),
            ", ",
            sQuote("alpha_v"),
            ") must be non-negative numeric."
        ),
        fixed = TRUE
    )

    # Prompt error when passing a matrix
    expect_error(
        moma_svd(X = X, lambda_v = matrix(1:12, 3)),
        paste0(
            "All penalty levels (",
            sQuote("lambda_u"),
            ", ",
            sQuote("lambda_v"),
            ", ",
            sQuote("alpha_u"),
            ", ",
            sQuote("alpha_v"),
            ") must be numeric."
        ),
        fixed = TRUE
    )

    expect_no_error(moma_svd(
        X = X,
        lambda_v = 1,
        lambda_u = 1
    ))
})


test_that("Arguments must be specified by name", {
    # lasso
    expect_error(
        lasso(FALSE),
        "Please specify the correct argument by name"
    )

    expect_no_error(lasso(non_negative = TRUE))

    # MCP
    expect_error(
        mcp(3),
        "Please specify the correct argument by name"
    )

    expect_error(
        mcp(3, non_negative = FALSE),
        "Please specify the correct argument by name"
    )
    expect_error(
        mcp(gamma = 3, TRUE),
        "Please specify the correct argument by name"
    )


    # grplasso
    expect_error(
        grplasso(factor(1)),
        "Please specify the correct argument by name"
    )
    expect_no_error(grplasso(g = factor(1)))


    # fused lasso
    expect_error(
        fusedlasso("path"),
        "Please specify the correct argument by name"
    )
    expect_no_error(fusedlasso(algo = "path"))

    # trend filtering
    expect_error(l1tf(1), "Please specify the correct argument by name")
    expect_no_error(l1tf(l1tf_k = 1))

    # sparse fused lasso
    expect_error(
        spfusedlasso(4),
        "Please specify the correct argument by name"
    )
    expect_no_error(spfusedlasso(lambda2 = 4))

    # clustering
    expect_error(cluster(diag(3)), "Please specify the correct argument by name")
    expect_no_error(cluster(w = diag(3)))

    # PG loop settings
    expect_error(
        moma_pg_settings(1e-10),
        "Please specify the correct argument by name"
    )
    expect_no_error(moma_pg_settings(EPS = 1e-10))
})
