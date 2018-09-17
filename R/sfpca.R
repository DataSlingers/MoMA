sfpca <- function(X,
                  # sparsity
                  P_v = "none",
                  P_u = "none",
                  lambda_v = 0,  # a vector or scalar, same for lambda_u, alpha_u/v
                  lambda_u = 0,
                  gamma_v=3,
                  gamma_u=3,
                  # non-negativity
                  nonneg_u = FALSE,
                  nonneg_v = FALSE,
                  # grouping
                  group_u = MOMA_EMPTYVEC,
                  group_v = MOMA_EMPTYVEC,
                  # sparse fused lasso
                  lambda2_u = 0,    # penalty on the abs value of parameters
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
                  # trend filtering
                  l1tf_k_u = 1,
                  l1tf_k_v = 1,
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
        moma_error("X must be a matrix.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]


    P_u <- toupper(P_u)
    P_v <- toupper(P_v)
    solver <- toupper(solver)

    alpha_u <- as.vector(alpha_u)
    alpha_v <- as.vector(alpha_v)
    lambda_u <- as.vector(lambda_u)
    lambda_v <- as.vector(lambda_v)

    Omega_u <- if(is.null(Omega_u)) diag(dim(X)[1]) else Omega_u
    Omega_v <- if(is.null(Omega_v)) diag(dim(X)[2]) else Omega_v

    return(cpp_sfpca(X = X,
                     w_v = w_v,w_u = w_u,
                     Omega_u = Omega_u,Omega_v = Omega_v,
                     alpha_u = alpha_u,alpha_v = alpha_v,
                     lambda_u = lambda_u,lambda_v = lambda_v,
                     P_u = P_u,P_v = P_v,gamma_u = gamma_u,gamma_v = gamma_v,
                     lambda2_u = lambda2_u, lambda2_v = lambda2_v,
                     ADMM_u = ADMM_u,ADMM_v = ADMM_v,
                     acc_u = acc_u,acc_v = acc_v,
                     prox_eps_u = prox_eps_u, prox_eps_v = prox_eps_v,
                     l1tf_k_u = l1tf_k_u, l1tf_k_v = l1tf_k_v,
                     nonneg_u = nonneg_u,nonneg_v = nonneg_v,
                     group_u = group_u,group_v = group_v,
                     EPS = EPS,MAX_ITER = MAX_ITER,
                     EPS_inner = EPS_inner,MAX_ITER_inner = MAX_ITER_inner,
                     solver = solver,
                     k = k))
}
