MOMA_EMPTYMAT <- matrix()
MOMA_EMPTYVEC <- vector(mode="numeric")

# This function checks
check_omega <- function(Omega,alpha,n){
    if(length(alpha) == 1 && alpha == 0){
        # discard the Omega matrix specified by users
        Omega <- diag(n)
    }
    else if(is.null(Omega)){
        # The user wants smooth penalty
        # but does not specify Omega matrix
        Omega <- second_diff_mat(n)
    }
    else{
        # Check validity of Omega if users speicify both alpha and Omega
        if(dim(Omega)[1] != dim(Omega)[2]){
            moma_error("Omega shoud be a square matrix: nrows = ",dim(Omega)[1],
                       ", ncols = ",dim(Omega)[2])
        }
        if(dim(Omega)[1] != n){
            moma_error("Omega shoud be a compatible matrix. It should be of ",
                       n,"x",n,
                       ", but is actually ",dim(Omega)[1],"x",dim(Omega)[1])
        }
    }
    return(Omega)
}

second_diff_mat <- function(n){
    return(crossprod(diff(diag(n))))
}

moma_svd <- function(
                    X,
                    u_sparsity=empty(),v_sparsity=empty(),lambda_u=0,lambda_v=0,    # lambda is a vector or scalar
                    Omega_u=NULL,Omega_v=NULL,alpha_u=0,alpha_v=0,                  # so is alpha
                    EPS = 1e-10, MAX_ITER = 1000,
                    EPS_inner = 1e-10,MAX_ITER_inner = 1e+5,
                    solver = "ista",
                    k = 1){

    # from scalar to vector
    alpha_u <- as.vector(alpha_u)
    alpha_v <- as.vector(alpha_v)
    lambda_u <- as.vector(lambda_u)
    lambda_v <- as.vector(lambda_v)

    df.arg.list <- list(
                        X = X,
                        P_v = "NONE",
                        P_u = "NONE",
                        lambda_u = lambda_u,
                        lambda_v = lambda_v,
                        gamma_u = 3,
                        gamma_v = 3,
                        # non-negativity
                        nonneg_u = FALSE,
                        nonneg_v = FALSE,
                        # grouping
                        group_u = MOMA_EMPTYVEC,
                        group_v = MOMA_EMPTYVEC,
                        lambda2_u = 0,
                        lambda2_v = 0,
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
                        alpha_u = alpha_u,
                        alpha_v = alpha_v,
                        # algorithm parameters
                        EPS = EPS,
                        MAX_ITER = MAX_ITER,
                        EPS_inner = EPS_inner,
                        MAX_ITER_inner = MAX_ITER_inner,
                        solver = toupper(solver),
                        k = k)
    if (!is.matrix(X)){
        moma_error("X must be a matrix.")
    }
    if (sum(!is.finite(X)) >= 1){
        moma_error("X must not have NaN, NA, or Inf.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]

    is_cv <- length(alpha_u) > 1 ||
              length(alpha_v) > 1 ||
              length(lambda_u) > 1 ||
              length(lambda_v) > 1
    # k must be 1 if alpha_u/v or lambda_u/v is of vector form
    if(is_cv && k != 1){
        moma_error("We don't support a range of parameters in finding a rank-k svd")
    }

    # Sparsity arguments
    if(!inherits(u_sparsity,"moma_sparsity") || !inherits(v_sparsity,"moma_sparsity")){
        moma_error("Sparse penalty should be of class '__moma_sp_'.
             Try using, for example, `u_sparsity = lasso()`.")
    }

    names(u_sparsity) <- if(length(u_sparsity) != 0) paste0(names(u_sparsity),"_u")
    # because paste0(NULL,"_u") =  "_u", check length(..) == 0 first
    names(v_sparsity) <- if(length(v_sparsity) != 0) paste0(names(v_sparsity),"_v")

    # We need u_sparsity to be a list
    arglist <- modifyList(df.arg.list,c(u_sparsity,v_sparsity))


    # Smoothness arguments
    arglist <- c(
                arglist,
                list(
                    Omega_u = check_omega(Omega_u,alpha_u,n),
                    Omega_v = check_omega(Omega_v,alpha_v,p))
                )

    if(is_cv){
        a <- do.call("cpp_sfpca_cv",arglist)
        class(a) <- "moma_svd_cv"
        return(a)
    }
    else{
        return(do.call("sfpca",arglist))
    }
}
