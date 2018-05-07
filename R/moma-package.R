#' MoMA: Modern Multivariate Analysis in R
#'
#' TODO
#'
#' See the package vignettes for details of the algorithm, as well as
#' comparisons to existing methods.
#'
#' @docType package
#' @name MoMA
#' @useDynLib MoMA
NULL

## Importing evalCpp is a hack to ensure that Rcpp is initialized before MoMA
## and may not be necessary if MoMA adds dependencies on other packages which
## will in turn initialize Rcpp
##
## Failing to have this import leads to a fairly nasty bug where
## MoMA:::.onAttach() fails and, more often than not, hangs indefinitely
#' @importFrom Rcpp evalCpp
NULL
