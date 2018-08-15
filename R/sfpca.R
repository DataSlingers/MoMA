MOMA_EMPTYMAT <- matrix()
MOMA_EMPTYVEC <- vector(mode="numeric")

sfpca <- function(X,
                  # sparsity
                  P_v = "lasso",
                  P_u = "lasso",
                  lambda_v = 0,
                  lambda_u = 0,
                  gamma=3,
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
                  Omega_u = NULL,
                  Omega_v = NULL,
                  alpha_u = 0,
                  alpha_v = 0,
                  # algorithm parameters
                  EPS = 1e-10,
                  MAX_ITER = 1000,
                  EPS_inner = 1e-10,
                  MAX_ITER_inner = 1e+5,
                  solver = "ista",
                  k = 1){
    if (!is.null(X) && !is.matrix(X)){
        stop("X must be a matrix.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]


    P_u <- toupper(P_u)
    P_v <- toupper(P_v)
    solver <- toupper(solver)
    # Smoothness arguments
    if(alpha_u == 0){
        # regardless of input Omega_u
        Omega_u <- diag(n)
    }
    else if(is.null(Omega_u)){
        # The user wants smooth penalty
        # but does not specify Omega matrix
        stop("here")
        Omega_u <- second.diff.mat(n)
    }
    if(alpha_v == 0){
        Omega_v <- diag(p)
    }
    else if(is.null(Omega_v)){
        stop("here")

        Omega_v <- second.diff.mat(p)
    }

    return(cpp_sfpca(X = X,
                     w_v = w_v,w_u = w_u,
                     Omega_u = Omega_u,Omega_v = Omega_v,
                     alpha_u = alpha_u,alpha_v = alpha_v,
                     lambda_u = lambda_u,lambda_v = lambda_v,
                     P_u = P_u,P_v = P_v,gamma = gamma,
                     ADMM_u = ADMM_u,ADMM_v = ADMM_v,
                     acc_u = acc_u,acc_v = acc_v,
                     prox_eps_u = prox_eps_u, prox_eps_v = prox_eps_v,
                     nonneg_u = nonneg_u,nonneg_v = nonneg_v,
                     group_u = group_u,group_v = group_v,
                     EPS = EPS,MAX_ITER = MAX_ITER,
                     EPS_inner = EPS_inner,MAX_ITER_inner = MAX_ITER_inner,
                     solver = solver,
                     k = k))
}
