sfpca <- function(X,
                  # sparsity
                  P_v = "none",
                  P_u = "none",
                  lambda_v = 0, # a vector or scalar, same for lambda_u, alpha_u/v
                  lambda_u = 0,
                  gamma_v = 3,
                  gamma_u = 3,
                  # non-negativity
                  nonneg_u = FALSE,
                  nonneg_v = FALSE,
                  # grouping
                  group_u = MOMA_EMPTYVEC,
                  group_v = MOMA_EMPTYVEC,
                  # sparse fused lasso
                  lambda2_u = 0, # penalty on the abs value of parameters
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
                  k = 1) {
    if (!is.null(X) && !is.matrix(X)) {
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

    Omega_u <- Omega_u %||% diag(dim(X)[1])
    Omega_v <- Omega_v %||% diag(dim(X)[2])

    prox_arg_list_u <- list(
        w = w_u,
        Omega = Omega_u,
        alpha = alpha_u,
        lambda = lambda_u,
        P = P_u,
        gamma = gamma_u,
        lambda2 = lambda2_u,
        ADMM = ADMM_u,
        acc = acc_u,
        prox_eps = prox_eps_u,
        l1tf_k = l1tf_k_u,
        nonneg = nonneg_u,
        group = group_u
    )
    prox_arg_list_v <- list(
        w = w_v,
        Omega = Omega_v,
        alpha = alpha_v,
        lambda = lambda_v,
        P = P_v,
        gamma = gamma_v,
        lambda2 = lambda2_v,
        ADMM = ADMM_v,
        acc = acc_v,
        prox_eps = prox_eps_v,
        l1tf_k = l1tf_k_v,
        nonneg = nonneg_v,
        group = group_v
    )
    return(cpp_moma_multi_rank(
        X = X,
        alpha_u = alpha_u, alpha_v = alpha_v,
        Omega_u = Omega_u, Omega_v = Omega_v,
        lambda_u = lambda_u, lambda_v = lambda_v,
        prox_arg_list_u = prox_arg_list_u,
        prox_arg_list_v = prox_arg_list_v,
        EPS = EPS, MAX_ITER = MAX_ITER,
        EPS_inner = EPS_inner, MAX_ITER_inner = MAX_ITER_inner,
        solver = solver,
        rank = k
    ))
}
