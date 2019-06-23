library(MoMA)
library(stringr)

expect_no_error <- function(object, ..., all = FALSE, info = NULL, label = NULL) {
    expect_error(object, regexp = NA, ..., all = all, info = info, label = label)
}

expect_no_warning <- function(object, ..., all = FALSE, info = NULL, label = NULL) {
    expect_warning(object, regexp = NA, ..., all = all, info = info, label = label)
}

expect_no_message <- function(object, ..., all = FALSE, info = NULL, label = NULL) {
    expect_message(object, regexp = NA, ..., all = all, info = info, label = label)
}

expect_str_contains <- function(object, expected, info = NULL, label = NULL) {
    if (!is.character(object)) object <- as.character(object)
    if (!is.character(expected)) expected <- as.character(expected)

    expect_true(all(str_detect(object, expected)),
        info = info, label = label
    )
}
