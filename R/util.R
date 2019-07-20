`%||%` <- function(lhs, rhs) {
    if (!is.null(lhs)) {
        lhs
    } else {
        rhs
    }
}

is_fin_numeric_scalar <- function(x) {
    (length(x) == 1L) && is.numeric(x) && (!is.na(x)) && is.finite(x)
}
is_valid_parameters <- function(x) {
    length(x) >= 1 && all(sapply(x, is_fin_numeric_scalar))
}
is_valid_select_str <- function(x) {
    is.character(x) && nchar(x) == 1 && x %in% c("b", "g")
}
# Check whether `x` is a boolean value
is_logical_scalar <- function(x) {
    return(is.logical(x) && (length(x) == 1) && !is.na(x))
}
