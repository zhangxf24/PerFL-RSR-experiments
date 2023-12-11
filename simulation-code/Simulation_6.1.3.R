### This simulation takes about 30 minutes for 
### one machine with 8 cores

rm(list = ls(all = TRUE))

library(hqreg)
library(Matrix)
library(doParallel)

# set path to the code folder
# path = setwd("..") 

source('simulation-code/Functions_evaluation.R')
source('simulation-code/Functions_regularizer.R')
source('simulation-code/Functions_other_models.R')
source('simulation-code/Functions_PerFL-RSR.R')

##### setting #####
M=10
np = c(50, 15)

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
iota <- 3 #mcp

r = (1+norm(ATA,"2"))*rho*tau + 1
H <- (r-rho*tau)*diag(1, nrow = nrow(ATA)) - rho*tau*ATA

##### tune lambda2 #####
simu <- function(lambda2, seed) {
  
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
  fit <- lapply(data, 
                FUN = function(XY) {
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  
                  ### Fitting by each unit with square loss ###
                  ## LASSO
                  fit.UHuberLasso <- hqreg(XY[,-(1:2)], XY[,1], lambda = lambdaSL)
                  BIC <- apply(fit.UHuberLasso$beta, 2, function(x) {
                    x[abs(x) < 1e-4] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberLasso <- fit.UHuberLasso$beta[,which.min(BIC)]
                  beta.UHuberLasso[abs(beta.UHuberLasso) < 1e-4] <- 0
                  
                  return(beta.UHuberLasso)
                })
  
  omega0 <- unlist(fit)
  
  fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                lambda1=0.00774, lambda2=lambda2, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                lambda1=1.668, lambda2=lambda2, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                lambda1=1.668, lambda2=lambda2, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  return(rbind(PMLasso = evallong(as.vector(fit1$omega),betam),
               PMMCP = evallong(as.vector(fit2$omega),betam),
               PMSCAD = evallong(as.vector(fit3$omega),betam)))
}

lambdaSL = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))

m=50

cl <- makeCluster(8)
registerDoParallel(cores=8)

for (t in 1:length(lambdaSL)) {
  results <- foreach(i = 1:m, .errorhandling = 'pass', 
                     .packages = c('hqreg')) %dopar% simu(lambdaSL[t], seed = 2023+i)
  
  save.image(paste0('simulation-results/lambda2_',t,'.Rdata'))
}      

stopCluster(cl)

##### tune lambda1 #####
simu <- function(lambda1, seed) {
  
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
  fit <- lapply(data, 
                FUN = function(XY) {
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  
                  ### Fitting by each unit with square loss ###
                  ## LASSO
                  fit.UHuberLasso <- hqreg(XY[,-(1:2)], XY[,1], lambda = lambdaSL)
                  BIC <- apply(fit.UHuberLasso$beta, 2, function(x) {
                    x[abs(x) < 1e-4] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberLasso <- fit.UHuberLasso$beta[,which.min(BIC)]
                  beta.UHuberLasso[abs(beta.UHuberLasso) < 1e-4] <- 0
                  
                  return(beta.UHuberLasso)
                })
  
  omega0 <- unlist(fit)
  
  fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                lambda1=lambda1, lambda2=0.001, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                lambda1=lambda1, lambda2=0.278, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                lambda1=lambda1, lambda2=0.278, 
                proportion = 1,
                omega0 = omega0, sig=rep(4,M), itermax = 200)
  
  return(rbind(PMLasso = evallong(as.vector(fit1$omega),betam),
               PMMCP = evallong(as.vector(fit2$omega),betam),
               PMSCAD = evallong(as.vector(fit3$omega),betam)))
}


lambdaSL = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))

m=50

cl <- makeCluster(8)
registerDoParallel(cores=8)

for (t in 1:length(lambdaSL)) {
  results <- foreach(i = 1:m, .errorhandling = 'pass',
                     .packages = c('hqreg')) %dopar% simu(lambdaSL[t], seed = 2023+i)
  
  save.image(paste0('simulation-results/lambda1_',t,'.Rdata'))

}      

stopCluster(cl)

##### tune sigma #####
simu <- function(sig, seed) {
  
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
  fit <- lapply(data, 
                FUN = function(XY) {
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  
                  ### Fitting by each unit with square loss ###
                  ## LASSO
                  fit.UHuberLasso <- hqreg(XY[,-(1:2)], XY[,1], lambda = lambdaSL)
                  BIC <- apply(fit.UHuberLasso$beta, 2, function(x) {
                    x[abs(x) < 1e-4] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberLasso <- fit.UHuberLasso$beta[,which.min(BIC)]
                  beta.UHuberLasso[abs(beta.UHuberLasso) < 1e-4] <- 0
                  
                  return(beta.UHuberLasso)
                })
  
  omega0 <- unlist(fit)
  
  fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                lambda1=0.00774, lambda2=0.001, 
                proportion = 1,
                omega0 = omega0, sig=rep(sig,M), itermax = 200)
  
  fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                lambda1=1.668, lambda2=0.278, 
                proportion = 1,
                omega0 = omega0, sig=rep(sig,M), itermax = 200)
  
  fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                lambda1=1.668, lambda2=0.278, 
                proportion = 1,
                omega0 = omega0, sig=rep(sig,M), itermax = 200)
  
  return(rbind(PMLasso = evallong(as.vector(fit1$omega),betam),
               PMMCP = evallong(as.vector(fit2$omega),betam),
               PMSCAD = evallong(as.vector(fit3$omega),betam)))
}

sigvs <- 1:7
m=50

cl <- makeCluster(8)
registerDoParallel(cores=8)

for (t in 1:length(sigvs)) {
  results <- foreach(i = 1:m, .errorhandling = 'pass', 
                     .packages = c('hqreg')) %dopar% simu(sigvs[t], seed = 2023+i)
  
  save.image(paste0('simulation-results/sigma_',sigvs[t],'.Rdata'))

}      

stopCluster(cl)