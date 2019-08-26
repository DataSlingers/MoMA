context("Tests for Util.R")

test_that("is_* functions", {
    expect_true(all(
        is_finite_numeric_scalar(1),
        !is_finite_numeric_scalar(c()),
        is_finite_numeric_scalar(c(1)),
        !is_finite_numeric_scalar("a"),
        !is_finite_numeric_scalar(Inf),
        !is_finite_numeric_scalar(NaN),
        !is_finite_numeric_scalar(NA)
    ))

    expect_true(all(
        # compatible with `is_finite_numeric_scalar`
        is_valid_parameters(1),
        !is_valid_parameters(c()),
        is_valid_parameters(c(1)),
        !is_valid_parameters("a"),
        !is_valid_parameters(Inf),
        !is_valid_parameters(NaN),
        !is_valid_parameters(NA),

        !is_valid_parameters(list()),
        !is_valid_parameters(c("a", 1)),
        !is_valid_parameters(list(1, list(1, 2, 3))),

        !is_valid_parameters(c(1, Inf)),
        !is_valid_parameters(c(Inf, 1)),
        !is_valid_parameters(c(1, NA)),
        !is_valid_parameters(c(NA)),
        !is_valid_parameters(c(1, NaN)),
        !is_valid_parameters(c(NaN, 1)),

        is_valid_parameters(seq(0, 10)),
        is_valid_parameters(matrix(c(1, 2, 3, 4), 2)),
        is_valid_parameters(c(1, 2, 3)),
        is_valid_parameters(list(1, 2, 3))
    ))


    expect_true(all(
        !is_valid_select_str(1),
        !is_valid_select_str("a"),
        is_valid_select_str("g"),
        is_valid_select_str("b"),
        !is_valid_select_str("bb"),
        !is_valid_select_str("")
    ))
})
