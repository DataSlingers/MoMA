`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs)) {
        lhs
    } else {
        rhs
    }
}

is_fin_numeric_scalar  <- function(x) {(length(x) == 1L) && is.numeric(x) && (!is.na(x)) && is.finite(x)}
is_valid_parameters <- function(x) {length(x) >= 1 && all(sapply(x, is_fin_numeric_scalar))}
