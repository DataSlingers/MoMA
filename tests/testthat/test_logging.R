library(stringr)
context("Logging Tests")

test_that("Logging controls work", {
    expect_error(moma_logger_level("BAD LEVEL"))

    moma_logger_level("INFO")
    expect_equal("INFO", moma_logger_level())

    moma_logger_level("MESSAGE")
})

test_that("INFO and DEBUG message print as expected", {
    moma_logger_level("MESSAGE")

    expect_silent(MoMA:::moma_info("A message"))
    expect_silent(MoMA:::moma_debug("A message"))

    moma_logger_level("DEBUG")

    expect_output(MoMA:::moma_info("A message"), "[INFO]")
    expect_output(MoMA:::moma_info("A message"), "A message")

    expect_output(MoMA:::moma_debug("The message"), "[DEBUG]")
    expect_output(MoMA:::moma_debug("The message"), "The message")

    moma_logger_level("MESSAGE")
})

test_that("Supressing messages works", {
    # At INFO level, everything is shown in R
    moma_logger_level("INFO")

    expect_error(MoMA:::moma_error("ERROR"))
    expect_warning(MoMA:::moma_warning("WARNING"))
    expect_message(MoMA:::moma_message("MESSAGE"))

    # At MESSAGE level, everything is shown in R
    moma_logger_level("MESSAGE")

    expect_error(MoMA:::moma_error("ERROR"))
    expect_warning(MoMA:::moma_warning("WARNING"))
    expect_message(MoMA:::moma_message("MESSAGE"))

    # At WARNING level, we don't get a message
    moma_logger_level("WARNING")

    expect_error(MoMA:::moma_error("ERROR"))
    expect_warning(MoMA:::moma_warning("WARNING"))
    expect_no_message(MoMA:::moma_message("MESSAGE"))

    # At ERROR level, we don't get a message or warning
    moma_logger_level("ERROR")

    expect_error(MoMA:::moma_error("ERROR"))
    expect_no_warning(MoMA:::moma_warning("WARNING"))
    expect_no_message(MoMA:::moma_message("MESSAGE"))

    moma_logger_level("MESSAGE")
})

test_that("No extra newlines", {
    moma_logger_level("DEBUG")

    e <- tryCatch(MoMA:::moma_error("MY ERROR"), error = identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(MoMA:::moma_warning("MY WARNING"), warning = identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(MoMA:::moma_message("MY MESSAGE"), message = identity)
    expect_equal(str_count(e$message, "\n"), 1)

    e <- tryCatch(MoMA:::moma_error("MY ERROR\nON TWO LINES"), error = identity)
    expect_equal(str_count(e$message, "\n"), 2)

    moma_logger_level("MESSAGE")
})

test_that("Function capture works at R level", {
    moma_logger_level("MESSAGE")

    f <- function(x) {
        MoMA:::moma_error("ERROR MESSAGE")
    }

    e <- tryCatch(f(), error = identity)
    expect_str_contains(e$message, "ERROR MESSAGE")
    expect_str_contains(e$message, "(Called from f)")
    expect_true(is.null(e$call))
    expect_true(is.null(e$cppstack))

    f <- function(x) {
        MoMA:::moma_error("ERROR MESSAGE", call = FALSE)
    }
    e <- tryCatch(f(), error = identity)

    expect_false(grepl("\\(Called from f\\)", e$message))

    f <- function(x) {
        MoMA:::moma_error("ERROR MESSAGE", call = "my func")
    }
    e <- tryCatch(f(), error = identity)

    expect_true(grepl("\\(Called from my func\\)", e$message))

    f <- function(x) {
        MoMA:::moma_warning("WARNING MESSAGE", call = FALSE)
    }
    e <- tryCatch(f(), warning = identity)

    expect_false(grepl("\\(Called from f\\)", e$message))

    f <- function(x) {
        MoMA:::moma_warning("WARNING MESSAGE", call = "my func")
    }
    e <- tryCatch(f(), warning = identity)

    expect_true(grepl("\\(Called from my func\\)", e$message))
})
