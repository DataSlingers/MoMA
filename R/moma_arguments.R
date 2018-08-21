empty <- function(){
    arglist <- list()
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

lasso <- function(nn = FALSE){
    if(!(nn %in% c(TRUE,FALSE))){
        stop("nn should be a boolean value.")
    }
    arglist <- list(nonneg = nn,P = "LASSO")
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

mcp <- function(gamma = 3, nn = FALSE){
    if(!(nn %in% c(TRUE,FALSE))){
        stop("nn should be a boolean value.")
    }
    arglist <- list(gamma = gamma, nonneg = nn, P = "MCP")
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

scad <- function(gamma = 3.7, nn = FALSE){
    if(!(nn %in% c(TRUE,FALSE))){
        stop("nn should be a boolean value.")
    }
    arglist <- list(gamma = gamma, nonneg = nn, P = "SCAD")
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

grplasso <- function(g, nn = FALSE){
    if(!(nn %in% c(TRUE,FALSE))){
        stop("nn should be a boolean value.")
    }
    if(class(g) != "factor"){
        stop("Please provide a factor for group lasso")
    }
    arglist <- list(group = g, P = "GRPLASSO", nonneg = nn)
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

fusedlasso <- function(){
    # fused lasso
    arglist <- list(P = "ORDEREDFUSED")
    class(arglist) <- "__moma_sp__"
    return(arglist)
}

cluster <- function(w = NULL,ADMM = FALSE,
                    acc = FALSE,
                    eps = 1e-10){
    # fused lasso
    if(!is.matrix(w) || is.null(w) || dim(w)[1] != dim(w)[2]){
        stop("`w` should be a square matrix.")
    }
    if(!isSymmetric(w)){
        warning("`w` is not symmetric. Only upper triangular half is used.")
    }
    arglist <- list(
                w = w, ADMM = ADMM,
                acc = acc, prox_eps = eps,
                P = "UNORDEREDFUSION")
    class(arglist) <- "__moma_sp__"
    return(arglist)
}
