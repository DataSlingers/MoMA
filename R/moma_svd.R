MOMA_EMPTYMAT <- matrix()
MOMA_EMPTYVEC <- vector(mode="numeric")

# This function checks
check_omega <- function(Omega,alpha,n){
    if(alpha == 0){
        # regardless of input Omeu
        Omega <- diag(n)
    }
    else if(is.null(Omega)){
        # The user wants smooth penalty
        # but does not specify Omega matrix
        Omega <- second.diff.mat(n)
    }else{
        if(dim(Omega)[1] != dim(Omega)[2] || dim(Omega)[1] != n){
            stop("Omega shoud be a compatible square matrix.")
        }
    }
    return(Omega)
}

second.diff.mat <- function(n){
    a <- diff(diag(n))
    return(t(a)%*%a)
}

moma_svd <- function(
                    X,
                    usp=empty(),vsp=empty(),lamu=0,lamv=0,
                    Omeu=NULL,Omev=NULL,alu=0,alv=0,
                    EPS = 1e-10, MAX_ITER = 1000,
                    EPS_inner = 1e-10,MAX_ITER_inner = 1e+5,
                    solver = "ista",
                    k = 1){

    df.arg.list <- list(
                        X = X,
                        P_v = "NONE",
                        P_u = "NONE",
                        lambda_v = lamu,
                        lambda_u = lamv,
                        gamma_u = 3,
                        gamma_v = 3,
                        # non-negativity
                        nonneg_u = FALSE,
                        nonneg_v = FALSE,
                        # grouping
                        group_u = MOMA_EMPTYVEC,
                        group_v = MOMA_EMPTYVEC,
                        # unordered fusion
                        w_u = MOMA_EMPTYMAT,
                        w_v = MOMA_EMPTYMAT,
                        ADMM_u = FALSE,
                        ADMM_v = FALSE,
                        acc_u = FALSE,
                        acc_v = FALSE,
                        prox_eps_u = 1e-10,
                        prox_eps_v = 1e-10,
                        # smoothness
                        alpha_u = alu,
                        alpha_v = alv,
                        # algorithm parameters
                        EPS = EPS,
                        MAX_ITER = MAX_ITER,
                        EPS_inner = EPS_inner,
                        MAX_ITER_inner = MAX_ITER_inner,
                        solver = solver,
                        k = k)
    if (!is.matrix(X)){
        stop("X must be a matrix.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]

    solver <- toupper(solver)

    # Sparsity arguments
    if(class(usp) != "__moma_sp__" || class(vsp) != "__moma_sp__"){
        stop("Sparse penalty should be of class '__moma_sp_'.
             Try using, for example, `usp = lasso()`.")
    }

    usp <- if(lamu == 0) list() else usp
    vsp <- if(lamv == 0) list() else vsp

    names(usp) <- if(length(usp) != 0) paste0(names(usp),"_u")
    names(vsp) <- if(length(vsp) != 0) paste0(names(vsp),"_v")

    arglist <- modifyList(df.arg.list,c(usp,vsp))


    # Smoothness arguments
    arglist <- c(
                arglist,
                list(
                    Omega_u = check_omega(Omeu,alu,n),
                    Omega_v = check_omega(Omev,alv,p))
                )

    return(do.call("sfpca",arglist))
}
