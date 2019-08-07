library(MoMA)
X <- rbind(iris3[, , 1], iris3[, , 2], iris3[, , 3])
Iris <- data.frame(X,
    Sp = rep(c("s", "c", "v"), rep(50, 3))
)
a <- SFPCA$new(
    X = Iris[, 1:4], rank = 2,
    v_sparsity = lasso(), lambda_v = seq(0, 0.5, 0.1),
    Omega_v = diag(4), alpha_v = seq(0, 0.5, 0.1)
)
a$plot()
