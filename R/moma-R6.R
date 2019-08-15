#' R6 objects for storing and accessing the results of SFPCA / SFLDA /
#' SFCCA
#'
#' During initialization of an \code{SFPCA} object, \code{R}
#' calls the \code{C++}-side function, \code{cpp_multirank_BIC_grid_search}, and
#' wraps the results returned. The \code{SFPCA} object also records penalty levels
#' and selection schemes of tuning parameters. Several helper
#' methods are provivded to facilitate access to results.
#' Initialization is delegated to \code{\link{moma_sfpca}}.
#' @seealso \code{\link{moma_sfpca}},
#' \code{\link{moma_sflda}},
#' \code{\link{moma_sfcca}}
#'
#' @section Members:
#'
#' \describe{
#'   \item{
#'      \code{center,scale}
#'   }{
#'      The attributes "\code{scaled:center}" and "\code{scaled:scale}" of function \code{scale}.
#'      The numeric centering and scalings used (if any) of the data matrix.
#'  }
#'
#'  \item{
#'      \code{grid_result}
#'  }{
#'      A 5-D list containing the results evaluated on the parameter grid.
#'  }
#'
#'  \item{
#'      \code{select_scheme_list}
#'  }{
#'      A list with elements 
#'      \code{select_scheme_alpha_u}, \code{select_scheme_alpha_v},
#'      \code{select_scheme_lambda_u}, \code{select_scheme_lambda_v}.
#'      Each of them is either 0 or 1. 0 stands for grid search
#'      and 1 stands for BIC search. Please see the \code{select_scheme}
#'      argument in the function \code{moma_sfpca}.
#'  }
#' }
#'
#' @section Methods:
#' \describe{
#'      \item{
#'          \code{get_mat_by_index}
#'      }{
#'         \describe{
#'              \item{
#'                  Arguments
#'              }{
#'                  \code{alpha_u}, \code{alpha_v}, \code{lambda_u}, \code{lambda_v}.
#'              }
#'              \item{}{
#'                  Indices of the parameters in the paramter grid, which have 
#'                  been specified during initialization.
#'              }
#'              \item{
#'                  Functionality
#'              }{
#'                  Obtain the right and left sigular penalized vectors, which are packed
#'                  into matrices \code{U} and \code{V}.
#'              }
#'          }
#'      }
#'      \item{
#'          \code{get_mat_by_index}
#'      }{
#'         \describe{
#'              \item{
#'                  Arguments
#'              }{
#'                  \code{alpha_u}, \code{alpha_v}, \code{lambda_u}, \code{lambda_v}.
#'              }
#'              \item{}{
#'                  Indices of the parameters in the paramter grid, which have 
#'                  been specified during initialization.
#'              }
#'              \item{
#'                  Functionality
#'              }{
#'                  Obtain the right and left sigular penalized vectors, which are packed
#'                  into matrices \code{U} and \code{V}.
#'              }
#'          }
#'      }
#'  }
#' @name moma_R6
NULL
