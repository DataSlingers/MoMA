#-------------------
# Util function
#-------------------
norm_vec <- function(x) sqrt(sum(x^2))

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

uni <- function(n){
    u_1 <- as.vector(rnorm(n))
    return(u_1/norm_vec(u_1))
}

#-------------------
# Generate data
#-------------------
set.seed(42)
n <- 199
p <- 200
ind <- as.vector(seq(p))
u_1 <- uni(n)
u_2 <- uni(n)
u_3 <- uni(n)
eps <- matrix(rnorm(n*p),n,p)
eps <- eps/20

# Sinusoidal
v_1 <- sin((ind+15)*pi/17);v_1[floor(7/20*p):p]=0;
v_1 <- v_1/norm_vec(v_1);
# Gaussian-modulated sinusoidal
v_2 <- as.vector(exp(-(ind-100)^2/650)*sin((ind-100)*2*pi/21));
v_2[0:floor(7/20*p)]=0;v_2[floor(130/200*p):p] = 0;
v_2 <- v_2/norm_vec(v_2);
# Sinusoidal
v_3 <- sin((ind-40)*pi/30);v_3[0:floor(130/200*p)]=0;
v_3 <- v_3/norm_vec(v_3);

plot(v_1,type = 'l',ylim=c(-0.3,0.3));
lines(v_2,col='blue');
lines(v_3,col='red')

X <-  n/4*u_1 %*% t(v_1) +eps + n/5*u_2 %*% t(v_2) #+ n/6 * u_3 %*% t(u_3)
# Noise proportion
print(norm(X) /norm(eps))
O_u <- SSD(n)
O_v <- SSD(p)
O_u <- diag(n)
#-------------------
# Demo 1, signal recovery
#-------------------
s <- svd(X)
plot(s$v[,1],type="l")
res1 <- sfpca(X,
              O_u,O_v,
              0,0,
              lambda_u=7,lambda_v=8,
              "SCAD","SCAD",
              gamma=3.7,
              EPS=1e-9,MAX_ITER = 1e+5,
              solver='ISTA')
plot(res1$v,type='l')
res2 <- sfpca(res1$DeflatedX,
              O_u,O_v,
              0,0,
              lambda_u=1,lambda_v=3,
              "LASSO","LASSO",
              gamma=3.7,
              EPS=1e-9,MAX_ITER=1e+5,
              solver='ISTA')
plot(res2$v,type='l')


res3 <- sfpca(res2$DeflatedX,
              O_u,O_v,
              alpha_u=1,alpha_v=1,

              lambda_u=6.5,lambda_v=1,
              "LASSO","LASSO",

              EPS=1e-9,MAX_ITER=1e+5,solver="ISTA")

plot(res3$v,type='l')


#-------------------
# Demo 2, effects of different penalty levels
#-------------------
sm_set = c(0.1,1,10,100)
sp_set = c(1,10,20,30)
par(mfrow=c(length(sm_set),length(sp_set)))
for(sm in sm_set){
    for(sp in sp_set){
        res <- sfpca(X,
                     O_u,O_v,
                     sm,sm,
                     lambda_u=sp,lambda_v=sp,
                     "LASSO","LASSO",
                     1e-9,1e+5)
        plot(res$v,type='l')
    }
}
title(main="Effects of penalty parameters",
      xlab="Sparsity", ylab="Smoothness")


W = matrix(c(0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,0),4,4)
I = diag(4)
C = I - 0.3*W
V = 4 * solve(t(C) %*% C)
image(-V)

print(xtable(V), floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)


## Generate animation
# par(mfrow = c(1,1))
# ani.options(interval=.05)
# sparse.p = c("l1","scad")
# cnt = 1
# for(l in seq(1,30,0.5)){
#     l1 <- sfpca("PCA",
#                 X=X,
#                 Y=matrix(runif(p*n),n,p),
#                 Omega_u=O_u,Omega_v=O_v,
#                 alpha_u=1,alpha_v=1,
#                 lambda_u=l,lambda_v=l,
#                 P_u="l1",P_v="l1",
#                 non_neg=0,
#                 scad_a=3.7,
#                 EPS=1e-9,
#                 MAX_ITER=1e+5,
#                 solver='ISTA',
#                 SVD = 1)
#     scad <- sfpca("PCA",
#                   X=X,
#                   Y=matrix(runif(p*n),n,p),
#                   Omega_u=O_u,Omega_v=O_v,
#                   alpha_u=1,alpha_v=1,
#                   lambda_u=l,lambda_v=l,
#                   P_u="scad",P_v="scad",
#                   non_neg=0,
#                   scad_a=3.7,
#                   EPS=1e-9,
#                   MAX_ITER=1e+5,
#                   solver='ISTA',
#                   SVD = 1)
#     l1nonneg <- sfpca("PCA",
#                       X=X,
#                       Y=matrix(runif(p*n),n,p),
#                       Omega_u=O_u,Omega_v=O_v,
#                       alpha_u=1,alpha_v=1,
#                       lambda_u=l,lambda_v=l,
#                       P_u="L1",P_v="L1",
#                       non_neg=1,
#                       scad_a=3.7,
#                       EPS=1e-9,
#                       MAX_ITER=1e+5,
#                       solver='ISTA',
#                       SVD = 1)
#     plot(l1$v,type='l',xlab=paste("scad//l1//Non-nega l1","lambda=",l),xlim=c(0,70))
#     lines(scad$v,col=cnt+2)
#     lines(l1nonneg$v,col=cnt+1)
#     save.image()
# }
# dev.off()
# system("convert  *.png example_1.gif")

