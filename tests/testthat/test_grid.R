context("`cpp_sfpca_grid`` is equivalent to run `cpp_sfpca` many times")

test_that("Using cpp_sfpca_grid is equivalent to run cpp_sfpca multiple times", {
    set.seed(332)
    n <- 7 # set n != p to avoid bugs
    p <- 11
    X = matrix(runif(n*p),n)

    # generate p.d. matrices
    O_v = crossprod(matrix(runif(p*p),p,p))
    O_u = crossprod(matrix(runif(n*n),n,n))

    # run tests
    # NOTE: there's no need to test for large
    # lambda's and alpha's because in those
    # cases u and v are zeros
    sp_set <- seq(0,3,0.5)
    sm_set <- seq(0,3,0.5)

    # WARNING: cannot add scad or mcp here
    # I guess because they are non-convex, so
    # there is slight difference in the results
    for(sptype in c(lasso,fusedlasso)){
        ista.cv <- moma_svd(X,
                            Omega_u=O_u,Omega_v=O_v,alpha_u=0,alpha_v=sm_set,
                            lambda_u=0,lambda_v=sp_set,u_sparsity=lasso(),v_sparsity=sptype(),
                            EPS=1e-14,MAX_ITER = 1e+5,solve="ISTA",EPS_inner = 1e-9)
        fista.cv <- moma_svd(X,
                            Omega_u=O_u,Omega_v=O_v,alpha_u=0,alpha_v=sm_set,
                            lambda_u=0,lambda_v=sp_set,u_sparsity=lasso(),v_sparsity=sptype(),
                            EPS=1e-14,MAX_ITER = 1e+5,solve="FISTA",EPS_inner = 1e-9)
        cnt = 1
        for(sp in sp_set){
            for(sm in sm_set){
                ista <- moma_svd(X,
                              Omega_u=O_u,Omega_v=O_v,alpha_u=0,alpha_v=sm,
                              lambda_u=0,lambda_v=sp,u_sparsity=lasso(),v_sparsity=sptype(),
                              EPS=1e-14,MAX_ITER = 1e+5,solve="ISTA",EPS_inner = 1e-9)
                fista <- moma_svd(X,
                               Omega_u=O_u,Omega_v=O_v,alpha_u=0,alpha_v=sm,
                               lambda_u=0,lambda_v=sp,u_sparsity=lasso(),v_sparsity=sptype(),
                               EPS=1e-14,MAX_ITER = 1e+5,solve="FISTA",EPS_inner = 1e-9)

                # Cannot use expect_equal here due to numerical error
                expect_lte(sum((ista$v[,1]-ista.cv$v[,cnt])^2),1e-7)
                expect_lte(sum((fista$v[,1]-fista.cv$v[,cnt])^2),1e-7)
                cnt = cnt + 1

            }
        }
    }
})
