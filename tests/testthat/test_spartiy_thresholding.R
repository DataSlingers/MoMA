context("Thresolding Tests")

test_that("Logging controls work", {
    expect_error(moma_logger_level("BAD LEVEL"))

    moma_logger_level("INFO")
    expect_equal("INFO", moma_logger_level())

    moma_logger_level("MESSAGE")
})
