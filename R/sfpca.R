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
                  # smoothness
                  Omega_u = NULL,
                  Omega_v = NULL,
                  alpha_u = 0,
                  alpha_v = 0,
                  # algorithm parameters
                  EPS = 1e-10,
                  MAX_ITER = 1000,
                  solver = "ista"){
    if (!is.null(X) && !is.matrix(X)){
        stop("X must be a matrix.")
    }
    n <- dim(X)[1]
    p <- dim(X)[2]
    P_u <- toupper(P_u)
    P_v <- toupper(P_v)
    solver <- toupper(solver)
    if(is.null(Omega_u)){
        Omega_u <- diag(n)
    }
    if(is.null(Omega_v)){
        Omega_v <- diag(p)
    }
    if(is.null(P_v)){
        stop("s")
    }
    return(cpp_sfpca(X = X,
                     w_v = w_v,w_u = w_u,
                     Omega_u = Omega_u,Omega_v = Omega_v,
                     alpha_u = alpha_u,alpha_v = alpha_v,
                     lambda_u = lambda_u,lambda_v = lambda_v,
                     P_u = P_u,P_v = P_v,gamma = gamma,
                     ADMM_u = ADMM_u,ADMM_v = ADMM_v,
                     acc_u = acc_u,acc_v = acc_v,
                     nonneg_u = nonneg_u,nonneg_v = nonneg_v,
                     group_u = group_u,group_v = group_v,
                     EPS = EPS,MAX_ITER = MAX_ITER,solver = solver))
}
