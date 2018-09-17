MOMA_EMPTYMAT <- matrix()
MOMA_EMPTYVEC <- vector(mode="numeric")
MOMA_DEFAULT_PROX <- list(
                        P = "NONE",
                        gamma = 3,
                        # non-negativity
                        nonneg = FALSE,
                        # grouping
                        group = MOMA_EMPTYVEC,
                        lambda2 = 0,
                        # unordered fusion
                        w = MOMA_EMPTYMAT,
                        ADMM = FALSE,
                        acc = FALSE,
                        prox_eps = 1e-10,
                        # trend filtering
                        l1tf_k = 1)

# This function checks the validity of Omega and alpha
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

    all_para <- c(alpha_u,alpha_v,lambda_u,lambda_v)
    if(sum(all_para < 0 || !is.finite(all_para)) > 0){
        moma_error("All penalty levels (",
                    sQuote("lambda_u"),", ",
                    sQuote("lambda_v"),", ",
                    sQuote("alpha_u"),", ",
                    sQuote("alpha_v"),
                    ") must be non-negative numeric.")
    }

    # from scalar to vector
    alpha_u <- as.vector(alpha_u)
    alpha_v <- as.vector(alpha_v)
    lambda_u <- as.vector(lambda_u)
    lambda_v <- as.vector(lambda_v)

    df_arg_list <- list(
                        X = X,
                        lambda_u = lambda_u,
                        lambda_v = lambda_v,
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
    df_prox_arg_list_u <- MOMA_DEFAULT_PROX
    df_prox_arg_list_v <- MOMA_DEFAULT_PROX

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
        moma_error("Sparse penalty should be of class ",
                    sQuote("moma_sparsity"),
                    ". Try using, for example, `u_sparsity = lasso()`.")
    }

    # Smoothness arguments
    df_arg_list <- c(
                df_arg_list,
                list(
                    Omega_u = check_omega(Omega_u,alpha_u,n),
                    Omega_v = check_omega(Omega_v,alpha_v,p),
                    prox_arg_list_u = modifyList(df_prox_arg_list_u,u_sparsity),
                    prox_arg_list_v = modifyList(df_prox_arg_list_v,v_sparsity)))

    if(is_cv){
        a <- do.call("cpp_sfpca_grid",df_arg_list)
        class(a) <- "moma_svd_grid"
        return(a)
    }
    else{
        return(do.call("cpp_sfpca",df_arg_list))
    }
}
