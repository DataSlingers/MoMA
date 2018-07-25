context("Unordered Fusion Lasso Tests")

# Special case 1:
# Large lambda: lambdas in each group are
# the means of everything in the group
test_that("Find means of everything when lambda is Large enough and the graph is fully connected", {
    set.seed(43)
    rep <- 10
    large.lambda <- 10000
    for(p in c(5,100)){
        y <- runif(p)
        w <- matrix(rep(0,p*p),p,byrow = T)
        for(i in 1:p-1){
            for(j in (i+1):p){
                w[i,j] = runif(1)+1
            }
        }
        for(i in 1:rep){
            res.AMA = test_prox_fusion(y,large.lambda,w,ADMM=FALSE,acc=FALSE)
            res.AMA.acc = test_prox_fusion(y,large.lambda,w,ADMM=FALSE,acc=TRUE)
            res.ADMM = test_prox_fusion(y,large.lambda,w,ADMM=TRUE,acc=FALSE)
            for(j in c(res.AMA,
                       res.AMA.acc,
                       res.ADMM)){
                expect_equal(j,mean(y))
            }
            print(i)
        }
    }
})

test_that("Find means of connected components when lambda is Large enough", {
    set.seed(33)
    rep <- 100
    large.lambda <- 1000
    for(p in c(9)){
        w <- matrix(rep(0,p*p),p,byrow = T);
        w[1,2] = runif(1) + 1
        w[1,3] = runif(1) + 1
        w[4,5] = runif(1) + 1
        w[4,6] = runif(1) + 1
        w[7,8] = runif(1) + 1
        w[7,9] = runif(1) + 1
        for(i in 1:rep){
            y <- 10 * runif(p)
            res.AMA.unacc = test_prox_fusion(y,large.lambda,w,ADMM=FALSE,acc=TRUE)
            res.AMA.acc = test_prox_fusion(y,large.lambda,w,ADMM=FALSE,acc=TRUE)
            res.ADMM = test_prox_fusion(y,large.lambda,w,ADMM=TRUE,acc=FALSE)
            for(res in list(res.AMA.unacc,
                            res.AMA.acc,
                            res.ADMM)){
                expect_equal(res[1],mean(y[1:3]))
                expect_equal(res[2],mean(y[1:3]))
                expect_equal(res[3],mean(y[1:3]))
                expect_equal(res[4],mean(y[4:6]))
                expect_equal(res[5],mean(y[4:6]))
                expect_equal(res[6],mean(y[4:6]))
                expect_equal(res[7],mean(y[7:9]))
                expect_equal(res[8],mean(y[7:9]))
                expect_equal(res[9],mean(y[7:9]))}
        }
    }
})

# Special case 2:
# For every lambda: when only w_ij = 1 (j = i+1)
# it becomes ordered fused lasso
test_that("Ordered fused lasso when w_ij = 1 all j = i+1", {
    set.seed(33)
    rep <- 100
    for(p in c(3,20,100)){
        w <- matrix(rep(0,p*p),p,byrow = T); for(i in 1:(p-1)) w[i,i+1] = 1
        for(i in 1:rep){

            print(i)
            y <- 10 * runif(p)
            err.AMA.unacc = norm(test_prox_orderedfusion(y,1)-test_prox_fusion(y,1,w,ADMM=FALSE,acc=TRUE))
            err.AMA.acc = norm(test_prox_orderedfusion(y,1)-test_prox_fusion(y,1,w,ADMM=FALSE,acc=TRUE))
            err.ADMM = norm(test_prox_orderedfusion(y,1)-test_prox_fusion(y,1,w,ADMM=TRUE,acc=FALSE))
            for(err in c(err.AMA.unacc,
                         err.AMA.acc,
                         err.ADMM)){
                expect_lte(err,1e-4)
            }
        }
    }
})

# Special case 3:
# For evry lambda: when all i,j w_ij = 1
# the solution path contains no split.
