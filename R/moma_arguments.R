# Check whether `x` is a boolean value
is_logical_scalar <- function(x){
    return(is.logical(x) && (length(x) == 1) && !is.na(x))
}

empty <- function(){
    arglist <- list()
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

lasso <- function(non_negative = FALSE){
    if(!is_logical_scalar(non_negative)){
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    arglist <- list(nonneg = non_negative,P = "LASSO")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

mcp <- function(gamma = 3, non_negative = FALSE){
    if(!is_logical_scalar(non_negative)){
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if(gamma < 1){
        moma_error("Non-convexity of MCP should be larger than 1.")
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "MCP")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

scad <- function(gamma = 3.7, non_negative = FALSE){
    if(!is_logical_scalar(non_negative)){
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if(gamma < 2){
        moma_error("Non-convexity of SCAD should be larger than 2.")
    }
    arglist <- list(gamma = gamma, nonneg = non_negative, P = "SCAD")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

grplasso <- function(g, non_negative = FALSE){
    if(!is_logical_scalar(non_negative)){
        moma_error(sQuote("non_negative"), " should be a boolean value.")
    }
    if(class(g) != "factor"){
        moma_error("Please provide a factor for group lasso")
    }
    arglist <- list(group = g, P = "GRPLASSO", nonneg = non_negative)
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

fusedlasso <- function(){
    # fused lasso
    arglist <- list(P = "ORDEREDFUSED")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}

cluster <- function(w = NULL,ADMM = FALSE,
                    acc = FALSE,
                    eps = 1e-10){
    # fused lasso
    if(!is.matrix(w) || is.null(w) || dim(w)[1] != dim(w)[2]){
        moma_error("`w` should be a square matrix.")
    }
    if(!isSymmetric(w)){
        moma_warning("`w` is not symmetric. Only upper triangular half is used.")
    }
    arglist <- list(
                w = w, ADMM = ADMM,
                acc = acc, prox_eps = eps,
                P = "UNORDEREDFUSION")
    class(arglist) <- "moma_sparsity"
    return(arglist)
}
