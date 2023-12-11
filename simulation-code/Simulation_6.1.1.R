### This simulation takes about 100 minutes for 
### one machine with 8 cores 15:50

library(ncvreg)
library(ILAMM)
library(Matrix)
library(doParallel)

# set path to the code folder
# path = setwd("..") 

source('simulation-code/Functions_evaluation.R')
source('simulation-code/Functions_regularizer.R')
source('simulation-code/Functions_other_models.R')
source('simulation-code/Functions_PerFL-RSR.R')

##### main simulation function #####

simu <- function(np, seed) {
  
  set.seed(seed)
  ##### data generation #####
  
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
  Ufit <- Unitfit(data)
  
  ##### FedAvg #####
  sigv = rep(4,M)
  fit2_avg <- fedavg(data, method = "Huber", regularizer = "MCP",
                     omega0 = Ufit$beta[,5], sig = sigv,
                     lambda = 0.1, proportion=1, itermax = 200)
  fit3_avg <- fedavg(data, method = "Huber", regularizer = "SCAD",
                     omega0 = Ufit$beta[,5], sig = sigv,
                     lambda = 0.1, proportion=1, itermax = 200)

  ##### Proposed method with fixed lambdas #####
  fit2_fix <- PMfit(data, method = "Huber", regularizer = "MCP",
                    lambda1=1.668, lambda2=0.278, proportion = 1,
                    omega0 = Ufit$beta[,5], sigv, itermax = 200)
  fit3_fix <- PMfit(data, method = "Huber", regularizer = "SCAD",
                    lambda1=1.668, lambda2=0.278, proportion = 1,
                    omega0 = Ufit$beta[,5], sigv, itermax = 200)
  
  ##### BIC #####
  
  sig <- max(apply(Ufit$tau.min, 1, min))
  sigv <- rep(sig,M)
  
  lambdaSL = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))
  lambd_ratio = c(3,6,9)
  nratio <- length(lambd_ratio)

  fit_MCP <- fit_SCAD <- list()
  BIC_MCP <- BIC_SCAD <- NULL
  for (i in 1:length(lambdaSL)) {
    for (j in 1:nratio) {
      fit2 <- PMfit(data, method = "Huber", regularizer = "MCP",
                    lambda1=lambd_ratio[j]*lambdaSL[i], 
                    lambda2=lambdaSL[i], proportion = 1,
                    omega0 = Ufit$beta[,5], sigv, itermax = 200)
      fit3 <- PMfit(data, method = "Huber", regularizer = "SCAD",
                    lambda1=lambd_ratio[j]*lambdaSL[i], 
                    lambda2=lambdaSL[i], proportion = 1,
                    omega0 = Ufit$beta[,5], sigv, itermax = 200)
      
      omega2 <- matrix(fit2$omega, ncol = M)
      omega3 <- matrix(fit3$omega, ncol = M)
      
      BIC2 <- BIC3 <- 0
      for (m in 1:M) {
        resid2 <- data[[m]][,1] - data[[m]][,-1] %*% omega2[,m]
        resid3 <- data[[m]][,1] - data[[m]][,-1] %*% omega3[,m]
        
        BIC2 <- BIC2 + 
          sum(ifelse(abs(resid2) <= sigv[m], 
                     resid2^2/2, sigv[m]*abs(resid2)-sigv[m]^2/2))/M/np[1] 
        BIC3 <- BIC3 + 
          sum(ifelse(abs(resid3) <= sigv[m], 
                     resid3^2/2, sigv[m]*abs(resid3)-sigv[m]^2/2))/M/np[1] 
      }
      
      BIC2 <- BIC2 + np[2]*log(np[1]*M)/np[1]/M*
        sum(apply(omega2, 1, function(x) length(unique(x))))
      BIC3 <- BIC3 + np[2]*log(np[1]*M)/np[1]/M*
        sum(apply(omega3, 1, function(x) length(unique(x))))
      
      BIC_MCP[(i-1)*nratio+j] <- BIC2
      BIC_SCAD[(i-1)*nratio+j] <- BIC3
        
      fit_MCP[[(i-1)*nratio+j]] <- fit2$omega
      fit_SCAD[[(i-1)*nratio+j]] <- fit3$omega
    }
  }
  
  beta.PMMCP <- fit_MCP[[which.min(BIC_MCP)]]
  beta.PMSCAD <- fit_SCAD[[which.min(BIC_SCAD)]]
  
  i = which.min(BIC_MCP) %/% nratio+1
  j = which.min(BIC_MCP) %% nratio

  return(list(eval = rbind(Ufit$eval,
                           AvgMCP = evallong(rep(fit2_avg$beta,M),betam),
                           AvgSCAD = evallong(rep(fit3_avg$beta,M),betam),
                           PMMCP = evallong(as.vector(fit2_fix$omega),betam),
                           PMSCAD = evallong(as.vector(fit3_fix$omega),betam),
                           BICMCP = evallong(as.vector(beta.PMMCP),betam),
                           BICSCAD = evallong(as.vector(beta.PMSCAD),betam)),
              par = c(sigma = sig, i = i, j = j)))
}


##### setting 1 #####
M = 10
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
iota <- 3 

r = (1+norm(ATA,"2"))*rho*tau + 1
H <- (r-rho*tau)*diag(1, nrow = nrow(ATA)) - rho*tau*ATA

Sigma <- matrix(0.3, nrow = np[2]-1, ncol = np[2]-1)
diag(Sigma) <- 1

m <- 50

cl <- makeCluster(8)
registerDoParallel(cores=8)

results <- foreach(i = 1:m, .errorhandling = 'pass', 
                   .packages = c('ncvreg','ILAMM')) %dopar% simu(np, seed = 2023+i)

stopCluster(cl)

save.image('simulation-results/BIC_M10_n50.Rdata')

### setting 2
M = 10
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

Sigma <- matrix(0.3, nrow = np[2]-1, ncol = np[2]-1)
diag(Sigma) <- 1

m <- 50

cl <- makeCluster(8)
registerDoParallel(cores=8)

results <- foreach(i = 1:m, .errorhandling = 'pass', 
                   .packages = c('ncvreg','ILAMM')) %dopar% simu(np, seed = 2023+i)

stopCluster(cl)

save.image('simulation-results/BIC_M10_n100.Rdata')

### setting 3
M = 50
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

Sigma <- matrix(0.3, nrow = np[2]-1, ncol = np[2]-1)
diag(Sigma) <- 1

m <- 50

cl <- makeCluster(8)
registerDoParallel(cores=8)

results <- foreach(i = 1:m, .errorhandling = 'pass', 
                   .packages = c('ncvreg','ILAMM')) %dopar% simu(np, seed = 2023+i)

stopCluster(cl)

save.image('simulation-results/BIC_M50_n100.Rdata')






