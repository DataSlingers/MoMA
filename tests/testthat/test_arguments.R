context("UX / Argument parsing tests")
test_that("Test for arguments names", {
    # This test makes sure the lasso(), mcp(), etc.,
    # return arguments with
    # correct names.

    correct_args <- names(formals(sfpca))

    # Collect all the arguments
    test_args <- list()

    for(fun in c(lasso,
                 scad,
                 mcp,
                 fusedlasso,
                 l1tf
                 )){
        test_args <- c(test_args, names(fun()))

    }

    args <- names(grplasso(g = factor(rep(1:10))))
    test_args <- c(test_args, args)

    args <- names(cluster(w = matrix(1)))
    test_args <- c(test_args, args)

    # Test
    for(arg in test_args){
        expect_true(paste0(arg,"_u") %in% correct_args)
        expect_true(paste0(arg,"_v") %in% correct_args)
    }
})

test_that("Prompt errors encountering inappropriate arguments", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))


    # Wrong non-convexity arguments
    expect_error(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=scad(1),lambda_u = 3),
                 paste0("Non-convexity parameter of SCAD (",sQuote("gamma"),") must be larger than 2."),fixed=TRUE)
    expect_error(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=mcp(0.9),lambda_u = 3),
                 paste0("Non-convexity parameter of MCP (",sQuote("gamma"),") must be larger than 1."),fixed=TRUE)


    # Wrong grouping dimension in group lasso
    expect_error(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=grplasso(factor(1)),lambda_u = 3),
                 "Wrong dimension: length(group) != dim(x).",fixed=TRUE)
    expect_error(grplasso(matrix(1)),
                 "Please provide a vector as an indicator of grouping. (Called from grplasso)",fixed=TRUE)


    # Wrong weight matrix dimension in cluster penalty
    expect_error(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=cluster(matrix(1)),lambda_u = 3),
                 "Wrong dimension: dim(weight matrix) != dim(x).",fixed=TRUE)


    # Omega has wrong dimension
    expect_error(moma_svd(matrix(runif(12),3,4),
                          Omega_u = matrix(c(1,2),1,2),alpha_u=2),
                 "Omega shoud be a square matrix: nrows = 1, ncols = 2 (Called from check_omega)",fixed=TRUE)
    expect_error(moma_svd(matrix(runif(12),3,4),
                          Omega_u = matrix(c(1),1),alpha_u=2),
                 "Omega shoud be a compatible matrix. It should be of 3x3, but is actually 1x1 (Called from check_omega)",fixed=TRUE)


    # Prompt errors when users require rank-k svd and cross validation
    expect_error(moma_svd(matrix(runif(12),3,4),lambda_u=c(1,2,3),k=2),
                 "We don't support a range of parameters in finding a rank-k svd (Called from moma_svd)",fixed=TRUE)
    expect_error(moma_svd(matrix(runif(12),3,4),
                          lambda_u=c(1,2,3),
                          lambda_v = seq(10),
                          alpha_u = seq(10)),
                 "We only allow changing two parameters.",fixed=TRUE)
})


test_that("Correct prox match", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))

    expect_output(moma_svd(matrix(runif(12),3,4)),
                  "Initializing null proximal operator object")


    # lasso
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=lasso(),lambda_u = 3),
                  "Initializing Lasso proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=lasso(TRUE),lambda_u = 3),
                  "Initializing non-negative Lasso proximal operator object")


    # scad
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=scad(),lambda_u = 3),
                  "Initializing SCAD proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=scad(non_negative=TRUE),lambda_u = 3),
                  "Initializing non-negative SCAD proximal operator object")


    # mcp
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=mcp(),lambda_u = 3),
                  "Initializing MCP proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=mcp(non_negative=TRUE),lambda_u = 3),
                  "Initializing non-negative MCP proximal operator object")


    # group
    expect_output(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=grplasso(factor(seq(3))),lambda_u = 3),
                 "Initializing group lasso proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=grplasso(factor(seq(3)),non_negative=TRUE),lambda_u = 3),
                  "Initializing non-negative group lasso proximal operator object")


    # L1 linear trend filtering
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=l1tf(),lambda_u = 3),
                  "Initializing a L1 linear trend filtering proximal operator object")
    # sparse fused lasso
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=spfusedlasso(lambda2=3),lambda_u = 3),
                  "Initializing a sparse fused lasso proximal operator object")
    # fused lasso
    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=fusedlasso(),lambda_u = 3),
                  "Initializing a ordered fusion lasso proximal operator object")

    # cluster penalty
    expect_output(moma_svd(matrix(runif(12),3,4),
                          u_sparsity=cluster(diag(3)),lambda_u = 3),
                 "Initializing a fusion lasso proximal operator object")
})

test_that("Correct algorithm match", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")


    expect_output(moma_svd(matrix(runif(12),3,4),solver = "ista"),
                  "Initializing a ISTA solver")
    expect_output(moma_svd(matrix(runif(12),3,4),solver = "fista"),
                  "Initializing a FISTA solver")
    expect_output(moma_svd(matrix(runif(12),3,4),solver = "onestepista"),
                  "Initializing an one-step ISTA solver")

    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=cluster(diag(3),ADMM=TRUE),lambda_u = 3),
                  "Running ADMM")

    expect_output(moma_svd(matrix(runif(12),3,4),
                           u_sparsity=cluster(diag(3)),lambda_u = 3),
                  "Running AMA")

    on.exit(MoMA::moma_logger_level(old_logger_level))
})

test_that("Data matrix must be complete", {
    old_logger_level <- MoMA::moma_logger_level()
    MoMA::moma_logger_level("DEBUG")
    on.exit(MoMA::moma_logger_level(old_logger_level))

    X <- matrix(runif(12),3,4)
    X[2,1] <- NA
    expect_error(moma_svd(X = X),
                  "X must not have NaN, NA, or Inf. (Called from moma_svd)",fixed=TRUE)

    X = matrix(runif(12),3,4)
    X[3,2] <- Inf
    expect_error(moma_svd(X = X),
                  "X must not have NaN, NA, or Inf. (Called from moma_svd)",fixed=TRUE)

    X = matrix(runif(12),3,4)
    X[1,4] <- NaN
    expect_error(moma_svd(X = X),
                  "X must not have NaN, NA, or Inf. (Called from moma_svd)",fixed=TRUE)
})
