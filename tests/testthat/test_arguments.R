context("Connection tests")
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
                 fusedlasso)){
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
    moma_set_logger_level_cpp(0)

    # Wrong non-convexity arguments
    expect_error(moma_svd(matrix(runif(12),3,4),
                          usp=scad(1),lamu = 3),
                 "Gamma for SCAD should be larger than 2!",fixed=TRUE)
    expect_error(moma_svd(matrix(runif(12),3,4),
                          usp=mcp(0.9),lamu = 3),
                 "Gamma for MCP should be larger than 1!",fixed=TRUE)

    # Wrong grouping dimension in group lasso
    expect_error(moma_svd(matrix(runif(12),3,4),
                          usp=grplasso(factor(1)),lamu = 3),
                 "Wrong dimension: dim(group) != dim(x).",fixed=TRUE)
    expect_error(grplasso(matrix(1)),
                 "Please provide a factor for group lasso",fixed=TRUE)


    # Wrong weight matrix dimension in cluster penalty
    expect_error(moma_svd(matrix(runif(12),3,4),
                          usp=cluster(matrix(1)),lamu = 3),
                 "Wrong dimension: dim(weight matrix) != dim(x).",fixed=TRUE)


    # Omega has wrong dimension
    expect_error(moma_svd(matrix(runif(12),3,4),
                          Omeu = matrix(c(1,2),1,2),alu=2),
                 "Omega shoud be a compatible square matrix.",fixed=TRUE)
    moma_set_logger_level_cpp(LEVELS[old_logger_level])
})


test_that("Correct prox match", {
    old_logger_level <- MoMA::moma_logger_level()
    moma_set_logger_level_cpp(0)


    expect_output(moma_svd(matrix(runif(12),3,4)),
                  "Initializing null proximal operator object")


    # lasso
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=lasso(),lamu = 3),
                  "Initializing Lasso proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=lasso(TRUE),lamu = 3),
                  "Initializing non-negative Lasso proximal operator object")


    # scad
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=scad(),lamu = 3),
                  "Initializing SCAD proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=scad(nn=TRUE),lamu = 3),
                  "Initializing non-negative SCAD proximal operator object")


    # mcp
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=mcp(),lamu = 3),
                  "Initializing MCP proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=mcp(nn=TRUE),lamu = 3),
                  "Initializing non-negative MCP proximal operator object")


    # group
    expect_output(moma_svd(matrix(runif(12),3,4),
                          usp=grplasso(factor(seq(3))),lamu = 3),
                 "Initializing group lasso proximal operator object")
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=grplasso(factor(seq(3)),nn=TRUE),lamu = 3),
                  "Initializing non-negative group lasso proximal operator object")


    # fused lasso
    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=fusedlasso(),lamu = 3),
                  "Initializing a ordered fusion lasso proximal operator object")

    # cluster penalty
    expect_output(moma_svd(matrix(runif(12),3,4),
                          usp=cluster(diag(3)),lamu = 3),
                 "Initializing a fusion lasso proximal operator object")

    moma_set_logger_level_cpp(LEVELS[old_logger_level])
})


test_that("Correct algorithm match", {
    old_logger_level <- MoMA::moma_logger_level()
    moma_set_logger_level_cpp(0)


    expect_output(moma_svd(matrix(runif(12),3,4),solver = "ista"),
                  "Initializing a ISTA solver")
    expect_output(moma_svd(matrix(runif(12),3,4),solver = "fista"),
                  "Initializing a FISTA solver")
    expect_output(moma_svd(matrix(runif(12),3,4),solver = "onestepista"),
                  "Initializing an one-step ISTA solver")

    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=cluster(diag(3),ADMM=TRUE),lamu = 3),
                  "Running ADMM")

    expect_output(moma_svd(matrix(runif(12),3,4),
                           usp=cluster(diag(3)),lamu = 3),
                  "Running AMA")

    moma_set_logger_level_cpp(LEVELS[old_logger_level])
})
