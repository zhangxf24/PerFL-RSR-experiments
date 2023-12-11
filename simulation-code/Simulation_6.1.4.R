### This simulation takes about an hour for 
### one machine with 8 cores

rm(list = ls(all = TRUE))

library(hqreg)
library(ncvreg)
library(Matrix)
library(doParallel)

# set path to the code folder
# path = setwd("..") 

path = setwd("/Users/xfzhang/Documents/Research/HuberReg/Submissions/JASA submission/Revision/PerFL-RSR-experiments")


source('simulation-code/Functions_evaluation.R')
source('simulation-code/Functions_regularizer.R')
source('simulation-code/Functions_other_models.R')
source('simulation-code/Functions_PerFL-RSR.R')

##### tune proportion #####
props <- (1:10)/10

M=10
np = c(100, 15)

beta1 <- c(-5,10,rep(0,np[2]-2-1))
beta2 <- c(5,10,rep(0,np[2]-2-1))
betam <- matrix(c(rep(beta1,M/2),rep(beta2,M/2)), byrow = FALSE, ncol = M)
betam <- rbind(rep(0,M),betam)

i <- rep(1:(M-1),(M-1):1)
j <- numeric()
for (m in 2:M) j <- c(j,m:M)
nM <- M*(M-1)/2
E <- sparseMatrix(i=c(1:nM,1:nM),j=c(i,j),
                  x=c(rep(1,nM),rep(-1,nM)))
OMEGA <- kronecker(E,diag(1,np[2]))
ATA <- t(OMEGA) %*% OMEGA

rho <- 2
tau <- 1
iota <- 3 

r = (1+norm(ATA,"2"))*rho*tau + 1
H <- (r-rho*tau)*diag(1, nrow = nrow(ATA)) - rho*tau*ATA


simu <- function(seed) {
  
  set.seed(seed)
  ##### data generation #####
  Sigma <- matrix(0.3, nrow = np[2]-1, ncol = np[2]-1)
  diag(Sigma) <- 1
  
  data <- list()
  for (i in 1:M) {
    X <- matrix(rnorm(np[1] * (np[2]-1)), nrow = np[1]) %*% Sigma
    X <- scale(X)
    if(i <= 5) {
      y <- X %*% beta1 + rt(np[1], df=1.5)
    } else {
      y <- X %*% beta2 + rt(np[1], df=3)
    }
    X <- cbind(rep(1,np[1]), X)
    data[[i]] <- cbind(y,X)
  }
  
  ##### Fitting by each unit #####
  Ufit <- Unitfit2(data)
  
  ### tune proportion ###
  sigv = rep(4,M)
  
  res <- list()
  for (i in 1:length(props)) {
    ##### Subgroup with l2 loss #####
    fit1_l2 <- PMfit_simu(data, method = "l2", regularizer = "lasso",
                          lambda1=0.00774, lambda2=0.00774, 
                          proportion = props[i],
                          omega0 = Ufit$beta[,1], sigv, itermax = 200/props[i])
    fit2_l2 <- PMfit_simu(data, method = "l2", regularizer = "MCP",
                          lambda1=1.668, lambda2=0.278, 
                          proportion = props[i],
                          omega0 = Ufit$beta[,1], sigv, itermax = 200/props[i])
    fit3_l2 <- PMfit_simu(data, method = "l2", regularizer = "SCAD",
                          lambda1=1.668, lambda2=0.278, 
                          proportion = props[i],
                          omega0 = Ufit$beta[,1], sigv, itermax = 200/props[i])
    
    ##### Proposed method #####
    fit1 <- PMfit_simu(data, method = "Huber", regularizer = "lasso",
                       lambda1=0.00774, lambda2=0.001, 
                       proportion = props[i],
                       omega0 = Ufit$beta[,4], sigv, itermax = 200/props[i])
    fit2 <- PMfit_simu(data, method = "Huber", regularizer = "MCP",
                       lambda1=1.668, lambda2=0.278, 
                       proportion = props[i],
                       omega0 = Ufit$beta[,4], sigv, itermax = 200/props[i])
    fit3 <- PMfit_simu(data, method = "Huber", regularizer = "SCAD",
                       lambda1=1.668, lambda2=0.278, 
                       proportion = props[i],
                       omega0 = Ufit$beta[,4], sigv, itermax = 200/props[i])
    
    res[[i]] <- rbind(SqLasso = fit1_l2$mse,
                      SqMCP = fit2_l2$mse,
                      SqSCAD = fit3_l2$mse,
                      PMLasso = fit1$mse,
                      PMMCP = fit2$mse,
                      PMSCAD = fit3$mse)
  }
  return(res)
}

m <- 50

cl <- makeCluster(8)
registerDoParallel(cores=8)

results <- foreach(i = 1:m, .errorhandling = 'pass',
                   .packages = c('ncvreg','hqreg')) %dopar% simu(seed = 2023+i)

stopCluster(cl)


save.image('simulation-results/proportion.Rdata')



