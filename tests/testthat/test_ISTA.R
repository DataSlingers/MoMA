#-------------------
# Util function
#-------------------
norm_vec <- function(x) sqrt(sum(x^2))

# generate second order difference matrix
SSD <- function(n){
    a <- 6*diag(n)
    for(i in 1:n){
        for(j in 1:n){
            if(abs(i-j) == 1) a[i,j] = -3;
            if(abs(i-j) == 2) a[i,j] = 1;
        }
    }
    return(a);
}

# a vector with norm 1
uni <- function(n){
    u_1 <- as.vector(rnorm(n))
    return(u_1/norm_vec(u_1))
}


#-------------------
# Generate data
#-------------------
set.seed(42)
n <- 199 # set n != p to test bugs
p <- 200
ind <- as.vector(seq(p)) # index from 1 to p
u_1 <- uni(n)
u_2 <- uni(n)
u_3 <- uni(n)
eps <- matrix(rnorm(n*p),n,p)
eps <- eps/20

# Sinusoidal signal
v_1 <- sin((ind+15)*pi/17);v_1[floor(7/20*p):p]=0;
v_1 <- v_1/norm_vec(v_1);
# Gaussian-modulated sinusoidal signal
v_2 <- as.vector(exp(-(ind-100)^2/650)*sin((ind-100)*2*pi/21));
v_2[0:floor(7/20*p)]=0;v_2[floor(130/200*p):p] = 0;
v_2 <- v_2/norm_vec(v_2);
# Sinusoidal signal
v_3 <- sin((ind-40)*pi/30);v_3[0:floor(130/200*p)]=0;
v_3 <- v_3/norm_vec(v_3);

plot(v_1,type = 'l',ylim=c(-0.3,0.3));
lines(v_2,col='blue');
lines(v_3,col='red')

X <-  n/4*u_1 %*% t(v_1) +eps + n/5*u_2 %*% t(v_2) #+ n/6 * u_3 %*% t(u_3)

# Noise proportion
print(norm(X) /norm(eps))

# Roughness penalty matrix
O_u <- SSD(n)
O_v <- SSD(p)
I_u <- diag(n)
I_v <- diag(p)

## Export data and run Genevera's matlab code for comparison
# writeMat("files.mat",v_1=v_1,v_2=v_2,v_3=v_3,u_1=u_1,u_2=u_2,u_3=u_3,X=X,O_u=O_u,O_v=O_v)

#-------------------
# test_that
#-------------------
test_that("Input validity test",{
})

test_that("Equal to SVD when no penalty", {
    sfpca <- sfpca(X,
                  O_u,O_v, 0,0,
                  lambda_u=0,lambda_v=0,"LASSO","LASSO",
                  gamma=3.7,EPS=1e-9,MAX_ITER = 1e+5,)
    svdd <- svd(X)
    expect_equal(sum((svdd$v[,1]-sfpca$v)^2),0)
    expect_equal(sum((svdd$u[,1]-sfpca$u)^2),0)
    expect_equal(svdd$d[1],sfpca$d);
    expect_error(moma_logger_level("BAD LEVEL"))
})

test_that("Equal to analytic solution when only roughness penalty n", {
    # TODO
})

test_that("Closed form solution when no sparsity",{
    # TODO
})

#-------------------
# Visual test
#-------------------
#-------------------
# Demo 1, signal recovery
#-------------------
svd <- svd(X)
plot(s$v[,1],type="l")
res1 <- sfpca(X,
              O_u,O_v,0,0,
              lambda_u=5,lambda_v=5,"LASSO","LASSO",
              solver='ISTA')
plot(res1$v,type='l')


#-------------------
# Demo 2, effects of different penalty levels
#-------------------
sm_set = c(0.1,1,10,100)
sp_set = c(1,3,5,6)
par(mfrow=c(length(sm_set),length(sp_set)))
for(sm in sm_set){
    for(sp in sp_set){
        res <- sfpca(X,
                     O_u,O_v,sm,sm,
                     lambda_u=sp,lambda_v=sp,"LASSO","LASSO",
                     EPS=1e-9,MAX_ITER=1e+5)
        plot(res$v,type='l')
    }
}
title(main="Effects of penalty parameters",
      xlab="Sparsity", ylab="Smoothness")


# Generate animation
# install.packages("animation")
# library(animation)
# par(mfrow = c(1,1))
# ani.options(interval=.05)
# sparse.p = c("l1","scad")
# cnt = 1
# for(l in seq(1,10,0.5)){
#     l1 <- sfpca(
#                 X=X,
#                 Omega_u=O_u,Omega_v=O_v,0,0,
#                 lambda_u=l,lambda_v=l,"LASSO","LASSO",gamma=3.7,
#                 EPS=1e-9,MAX_ITER=1e+5, solver='ISTA')
#     scad <- sfpca(
#                 X=X,
#                 Omega_u=O_u,Omega_v=O_v,0,0,
#                 lambda_u=l,lambda_v=l,"SCAD","SCAD",gamma=3.7,
#                 EPS=1e-9,MAX_ITER=1e+5, solver='ISTA')
#     l1nonneg <- sfpca(
#                 X=X,
#                 Omega_u=O_u,Omega_v=O_v,0,0,
#                 lambda_u=l,lambda_v=l,"MCP","MCP",gamma=3.7,
#                 EPS=1e-9,MAX_ITER=1e+5, solver='ISTA')
#     plot(l1$v,type='l',xlab=paste("scad/l1/Non-nega l1","lambda=",l),xlim=c(0,70))
#     lines(scad$v,col=cnt+2)
#     lines(l1nonneg$v,col=cnt+1)
#     save.image()
# }
# dev.off()


