rm(list = ls(all = TRUE))

library(expm)
library(ncvreg)
library(ILAMM)
library(Matrix)

path = setwd("/Users/xfzhang/Dropbox/Research/HuberReg/Code")
path = setwd("/Users/xiaofei/Dropbox/Research/HuberReg/Code")
source("Functions.R")

##### data generation #####
M=10
np = c(20, 100)

Sigma <- matrix(0.3, nrow = np[2]-1, ncol = np[2]-1)
diag(Sigma) <- 1
beta1 <- c(-2,4,rep(0,np[2]-2-1))
beta2 <- c(2,4,rep(0,np[2]-2-1))
betam <- matrix(c(rep(beta1,M/2),rep(beta2,M/2)), byrow = FALSE, ncol = M)
betam <- rbind(rep(0,10),betam)

data <- list()
for (i in 1:M) {
  X <- matrix(rnorm(np[1] * (np[2]-1)), nrow = np[1]) %*% Sigma
  X <- scale(X)
  if(i <= 5) {
    y <- X %*% beta1 
    # + rnorm(np[1], 0.1)
    y <- y - mean(y)
  } else {
    y <- X %*% beta2 
    # + rnorm(np[1], 0.1)
    y <- y - mean(y)
  }
  X <- cbind(rep(1,np[1]), X)
  data[[i]] <- cbind(y,X)
}

##### Fitting #####

fit <- lapply(data, 
              FUN = function(XY) {
                lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                  exp(seq(log(10^1), log(10^-6), length = 20))
                
                ### Fitting by each unit with square loss ###
                
                ### Lasso
                fit.USqLasso <- ncvreg(XY[,-(1:2)], XY[,1],
                                       penalty = "lasso", lambda = lambdaSL)
                BIC <- apply(fit.USqLasso$beta, 2, function(x) {
                  x[abs(x) < 1e-4] <- 0
                  log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                    sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                })
                beta.USqLasso <- fit.USqLasso$beta[,which.min(BIC)]
                beta.USqLasso[abs(beta.USqLasso) < 1e-4] <- 0
                ### MCP
                fit.USqMCP<- ncvreg(XY[,-(1:2)], XY[,1],
                                       penalty = "MCP", lambda = lambdaSL)
                BIC <- apply(fit.USqMCP$beta, 2, function(x) {
                  x[abs(x) < 1e-4] <- 0
                  log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                    sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                })
                beta.USqMCP <- fit.USqMCP$beta[,which.min(BIC)]
                beta.USqMCP[abs(beta.USqMCP) < 1e-4] <- 0
                ### SCAD
                fit.USqSCAD <- ncvreg(XY[,-(1:2)], XY[,1],
                                       penalty = "SCAD", lambda = lambdaSL)
                BIC <- apply(fit.USqSCAD$beta, 2, function(x) {
                  x[abs(x) < 1e-4] <- 0
                  log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                    sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                })
                beta.USqSCAD <- fit.USqSCAD$beta[,which.min(BIC)]
                beta.USqSCAD[abs(beta.USqSCAD) < 1e-4] <- 0
                
                ### Fitting by each unit with Huber loss ###
                
                ### Lasso
                fit.UHuberLasso <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "Lasso", 
                                                  lSeq = lambdaSL, nfolds = 2)
                fit.UHuberMCP <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "MCP", 
                                                  lSeq = lambdaSL, nfolds = 2)
                fit.UHuberSCAD <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "SCAD", 
                                                  lSeq = lambdaSL, nfolds = 2)
                # plot(beta.USqLasso,fit.UHuberMCP$beta)
                # c(fit.UHuberLasso$tauMin, fit.UHuberMCP$tauMin, fit.UHuberSCAD$tauMin)
                
                return(list(USqLasso = beta.USqLasso,
                            USqMCP = beta.USqMCP,
                            USqSCAD = beta.USqSCAD,
                            UHuberLasso = fit.UHuberLasso$beta,
                            UHuberMCP = fit.UHuberMCP$beta,
                            UHuberSCAD = fit.UHuberSCAD$beta,
                            tau = c(fit.UHuberLasso$tauMin, fit.UHuberMCP$tauMin, fit.UHuberSCAD$tauMin)))
              })

beta.USqLasso <- c(sapply(fit,"[[",1))
beta.USqMCP <- c(sapply(fit,"[[",2))
beta.USqSCAD <- c(sapply(fit,"[[",3))
beta.UHuberLasso <- c(sapply(fit,"[[",4))
beta.UHuberMCP <- c(sapply(fit,"[[",5))
beta.UHuberSCAD <- c(sapply(fit,"[[",6))


eval.USqLasso <- evallong(beta.USqLasso, betam)
eval.USqMCP <- evallong(beta.USqMCP, betam)
eval.USqSCAD <- evallong(beta.USqSCAD, betam)
eval.UHuberLasso <- evallong(beta.UHuberLasso, betam)
eval.UHuberMCP <- evallong(beta.UHuberMCP, betam)
eval.UHuberSCAD <- evallong(beta.UHuberSCAD, betam)

rbind(eval.USqLasso, eval.USqMCP, eval.USqSCAD,
      eval.UHuberLasso, eval.UHuberMCP, eval.UHuberSCAD)

### Proposed method

i <- rep(1:(M-1),(M-1):1)
j <- numeric()
for (m in 2:M) j <- c(j,m:M)
nM <- M*(M-1)/2
E <- sparseMatrix(i=c(1:nM,1:nM),j=c(i,j),
                  x=c(rep(1,nM),rep(-1,nM)))
OMEGA <- kronecker(E,diag(1,np[2]))
A <- rbind(OMEGA, diag(1,M*np[2]))
rho <- (1+sqrt(33))/2
ATA <- t(A) %*% A
# range(eigen(ATA)$values)
# range(eigen(H)$values)
eta <- 0.1
r = norm(ATA,"2")*rho*eta+1/M
H <- r*diag(1, nrow = nrow(ATA)) - rho*eta*ATA
nu <- 3 #mcp
# nu <- 3.7 #scad

data <- data_simu

itermax = 2000
penalty <- "lasso"
lambda1 <- 0.1
lambda2 <- 0.001
lambda <- c(rep(lambda1,nrow(A) - M*np[2]),rep(lambda2,M*np[2]))
tau = 5

omega_usq <- UnitSq(data_simu, "lasso", betam)$omega
omega0 <- omega_usq
# omega <- rnorm(M*np[2],1)
delta <- A %*% omega0
gamma <- rep(0,nrow(A))
mse <- lagrg <- huberloss <- numeric()
diff <- iter <- 1
while (iter <= itermax & diff > 0.0001) {
  mse[iter] <- norm(omega0-c(betam), type = "2")
  huloss <- 0
  for (m in 1:M) {
    yy <- data[[m]][,1]
    XX <- data[[m]][,-1]
    beta <- omega0[(1:np[2])+(m-1)*np[2]]
    resid <- yy - XX %*% beta
    huloss <- huloss+sum(ifelse(abs(resid) <= tau, resid^2/2, tau*abs(resid)-tau^2/2))
  }
  huberloss[iter] <- huloss/M/np[1]
  lagrg[iter] <- huloss/M/np[1] + sum(abs(lambda*delta)) +
    sum(gamma*(A%*%omega0 - delta)) + rho/2*sum((A%*%omega0-delta)^2)
  
  td <- A %*% omega0 + gamma/rho
  if(penalty == 'MCP') {
    delta <- ifelse(abs(td) <= lambda*nu, 
                    ifelse(abs(td) <= lambda/rho, 0, sign(td)*(abs(td)-lambda/rho))/(1-1/nu/rho),td)
  } else if(penalty == "SCAD") {
    delta <- sapply(td, function(aa) {
      if(aa <= lambda+lambda/rho) {
        delta <- ifelse(aa <= lambda/rho, 0, sign(aa)*(abs(aa)-lambda/rho))
      } else if(aa <= nu*lambda) {
        delta <- ifelse(aa <= nu*lambda/rho/(nu-1), 0, sign(aa)*(abs(aa)-nu*lambda/rho/(nu-1)))/(1-1/(nu-1)/rho)
      } else {
        delta <- aa
      }})
  } else if(penalty == "lasso" | penalty == "Lasso") {
    delta <- ifelse(abs(td) <= lambda/rho, 0, sign(td)*(abs(td)-lambda/rho))
  }
  
  grad <- numeric()
  for (m in 1:M) {
    yy <- data[[m]][,1]
    XX <- data[[m]][,-1]
    beta <- omega0[(1:np[2])+(m-1)*np[2]]
    resid <- yy - XX %*% beta
    grad[(1:np[2])+(m-1)*np[2]] <- -t(ifelse(abs(resid) <= tau, resid, sign(resid)*tau)) %*% XX
  }
  
  omega <- H %*% omega0 / r -
    eta / r * (grad/M/np[1] - t(A) %*% (rho * delta - gamma))
  # 
  # omega[omega < 1e-6] <- 0
  # delta[delta < 1e-6] <- 0
  
  gamma <- gamma + rho*(A %*% omega - delta)
  
  diff <- norm(omega - omega0, "2")
  iter <- iter + 1
  omega0 <- omega
}
plot(mse, type = 'l')
plot(lagrg, type = 'l')
plot(huberloss, type = 'l')

omega <- as.vector(omega)
omega[abs(omega) < 1e-4] <- 0
delta[abs(delta) < 1e-4] <- 0
plot(omega_usq,omega)
# unique(omega)
length(unique(omega))
length(unique(delta))
rho/2*sum((A%*%omega-delta)^2)
omega[1:10]
length(unique(as.vector(A %*% omega_usq)))















