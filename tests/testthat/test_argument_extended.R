context("Test extended argument helpers")

test_that("moma_* works", {
    expect_true(all(
        moma_lasso()$sparsity_type$nonneg == FALSE,
        moma_lasso()$sparsity_type$P == "LASSO",
        moma_lasso()$sparsity_type$nonneg == FALSE,
        class(moma_lasso()) == "moma_sparsity_type",
        moma_lasso(non_negative = TRUE)$sparsity_type$nonneg == TRUE,
        all(c("lambda", "select_scheme", "non_negative") %in% names(formals(moma_lasso))),
        moma_lasso()$lambda == 0,
        moma_lasso()$select_scheme == "grid"
    ))

    expect_warning(
        moma_lasso(TRUE),
        "extra argument  will be disregarded"
    )


    expect_true(all(
        moma_spfusedlasso(lambda2 = 2, lambda = seq(0, 10))$sparsity_type$lambda2 == 2,
        moma_spfusedlasso(lambda2 = 2, lambda = seq(0, 10))$lambda == seq(0, 10)
    ))
})
