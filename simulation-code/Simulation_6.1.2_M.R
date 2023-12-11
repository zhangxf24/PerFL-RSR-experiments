### This simulation takes about 70 minutes for 
### one machine with 8 cores

rm(list = ls(all = TRUE))

library(hqreg)
library(ncvreg)
library(Matrix)
library(doParallel)

# set path to the code folder
# path = setwd("..") 

source('simulation-code/Functions_evaluation.R')
source('simulation-code/Functions_regularizer.R')
source('simulation-code/Functions_other_models.R')
source('simulation-code/Functions_PerFL-RSR.R')

##### tune t #####
M_list <- c(20,50,100)

for (pp in 1:length(M_list)) {
  M=M_list[pp]
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
      if(i <= M/2) {
        y <- X %*% beta1 + rt(np[1], df=1.5)
      } else {
        y <- X %*% beta2 + rt(np[1], df=3)
      }
      X <- cbind(rep(1,np[1]), X)
      data[[i]] <- cbind(y,X)
    }
    
    
    ##### Fitting by each unit #####
    Ufit <- Unitfit2(data)
    
    ##### Subgroup with l2 loss #####
    
    fit1_l2 <- PMfit(data, method = "l2", regularizer = "lasso",
                     lambda1=0.00774, lambda2=0.00774, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit2_l2 <- PMfit(data, method = "l2", regularizer = "MCP",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit3_l2 <- PMfit(data, method = "l2", regularizer = "SCAD",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    
    ##### Proposed method #####
    
    fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                  lambda1=0.00774, lambda2=0.001, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    
    return(rbind(Ufit$eval,
                 SqLasso = evallong(as.vector(fit1_l2$omega),betam),
                 SqMCP = evallong(as.vector(fit2_l2$omega),betam),
                 SqSCAD = evallong(as.vector(fit3_l2$omega),betam),
                 PMLasso = evallong(as.vector(fit1$omega),betam),
                 PMMCP = evallong(as.vector(fit2$omega),betam),
                 PMSCAD = evallong(as.vector(fit3$omega),betam)))
  }
  
  ### run ###
  
  m <- 50
  set.seed(2023)
  
  cl <- makeCluster(8)
  registerDoParallel(cores=8)
  
  results <- foreach(i = 1:m, .errorhandling = 'pass',
                     .packages = c('ncvreg','hqreg')) %dopar% simu(seed = 2023+i)
  
  stopCluster(cl)
  
  save.image(paste0('simulation-results/t_M',M_list[pp],'_n100_p15.Rdata'))

}

##### normal distribution #####
M_list <- c(20,50,100)

for (pp in 1:length(M_list)) {
  M=M_list[pp]
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
      if(i <= M/2) {
        y <- X %*% beta1 + rnorm(np[1], sd=1)
      } else {
        y <- X %*% beta2 + rnorm(np[1], sd=0.5)
      }
      X <- cbind(rep(1,np[1]), X)
      data[[i]] <- cbind(y,X)
    }
    
    
    ##### Fitting by each unit #####
    Ufit <- Unitfit2(data)
    
    ##### Subgroup with l2 loss #####
    
    fit1_l2 <- PMfit(data, method = "l2", regularizer = "lasso",
                     lambda1=0.00774, lambda2=0.00774, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit2_l2 <- PMfit(data, method = "l2", regularizer = "MCP",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit3_l2 <- PMfit(data, method = "l2", regularizer = "SCAD",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    
    ##### Proposed method #####
    
    fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                  lambda1=0.00774, lambda2=0.001, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    
    return(rbind(Ufit$eval,
                 SqLasso = evallong(as.vector(fit1_l2$omega),betam),
                 SqMCP = evallong(as.vector(fit2_l2$omega),betam),
                 SqSCAD = evallong(as.vector(fit3_l2$omega),betam),
                 PMLasso = evallong(as.vector(fit1$omega),betam),
                 PMMCP = evallong(as.vector(fit2$omega),betam),
                 PMSCAD = evallong(as.vector(fit3$omega),betam)))
  }
  
  ### run ###
  
  m <- 50
  
  cl <- makeCluster(8)
  registerDoParallel(cores=8)
  
  results <- foreach(i = 1:m, .errorhandling = 'pass',
                     .packages = c('ncvreg','hqreg')) %dopar% simu(seed = 2023+i)
  
  stopCluster(cl)
  
  save.image(paste0('simulation-results/normal_M',M_list[pp],'_n100_p15.Rdata'))
  
}

##### cauchy distribution #####
M_list <- c(20,50,100)

for (pp in 1:length(M_list)) {
  M=M_list[pp]
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
      if(i <= M/2) {
        y <- X %*% beta1 + rcauchy(np[1], scale = 1)
      } else {
        y <- X %*% beta2 + rcauchy(np[1], scale = 1.5)
      }
      X <- cbind(rep(1,np[1]), X)
      data[[i]] <- cbind(y,X)
    }
    
    
    ##### Fitting by each unit #####
    Ufit <- Unitfit2(data)
    
    ##### Subgroup with l2 loss #####
    
    fit1_l2 <- PMfit(data, method = "l2", regularizer = "lasso",
                     lambda1=0.00774, lambda2=0.00774, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit2_l2 <- PMfit(data, method = "l2", regularizer = "MCP",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    fit3_l2 <- PMfit(data, method = "l2", regularizer = "SCAD",
                     lambda1=1.668, lambda2=0.278, proportion = 1,
                     omega0 = Ufit$beta[,1], rep(4,M), itermax = 200)
    
    ##### Proposed method #####
    
    fit1 <- PMfit(data, method = "Huber", regularizer = "lasso",
                  lambda1=0.00774, lambda2=0.001, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                  lambda1=1.668, lambda2=0.278, proportion = 1,
                  omega0 = Ufit$beta[,4], rep(4,M), itermax = 200)
    
    return(rbind(Ufit$eval,
                 SqLasso = evallong(as.vector(fit1_l2$omega),betam),
                 SqMCP = evallong(as.vector(fit2_l2$omega),betam),
                 SqSCAD = evallong(as.vector(fit3_l2$omega),betam),
                 PMLasso = evallong(as.vector(fit1$omega),betam),
                 PMMCP = evallong(as.vector(fit2$omega),betam),
                 PMSCAD = evallong(as.vector(fit3$omega),betam)))
  }
  
  ### run ###
  
  m <- 50
  
  cl <- makeCluster(8)
  registerDoParallel(cores=8)
  
  results <- foreach(i = 1:m, .errorhandling = 'pass',
                     .packages = c('ncvreg','hqreg')) %dopar% simu(seed = 2023+i)
  
  stopCluster(cl)
  
  save.image(paste0('simulation-results/cauchy_M',M_list[pp],'_n100_p15.Rdata'))
  
}
