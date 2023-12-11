##### Fitting by each unit #####

Unitfit <- function(data) {
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
                  
                  ### MCP
                  fit.UHuberMCP <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "MCP", 
                                                  lSeq = lambdaSL, nfolds = 2)
                  
                  ### SCAD
                  fit.UHuberSCAD <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "SCAD", 
                                                   lSeq = lambdaSL, nfolds = 2)
                  
                  return(list(USqLasso = beta.USqLasso,
                              USqMCP = beta.USqMCP,
                              USqSCAD = beta.USqSCAD,
                              UHuberLasso =fit.UHuberLasso$beta,
                              UHuberMCP = fit.UHuberMCP$beta,
                              UHuberSCAD = fit.UHuberSCAD$beta,
                              tau.min = c(fit.UHuberMCP$tauMin, fit.UHuberSCAD$tauMin),
                              lambda.min = c(fit.UHuberMCP$lambdaMin, 
                                             fit.UHuberSCAD$lambdaMin)))
                })
  
  beta.Ufit <- NULL
  for (i in 1:6) {
    beta.Ufit <- cbind(beta.Ufit,c(sapply(fit,"[[",i)))  
  }
  eval.Ufit <- t(apply(beta.Ufit, 2, function(x) evallong(x, betam)))
  row.names(eval.Ufit) <- c('USqLasso','USqMCP','USqSCAD',
                            'UHuberLasso','UHuberMCP','UHuberSCAD')
  
  return(list(beta = beta.Ufit,
              eval = eval.Ufit, 
              tau.min = matrix(c(sapply(fit,"[[",7)), nrow = M, byrow = TRUE), 
              lambda.min = matrix(c(sapply(fit,"[[",8)), nrow = M, byrow = TRUE)))
}

### use BIC to choose sigma
Unitfit1 <- function(data,sig) {
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
                  fit.UHuberLasso <- matrix(NA, nrow = ncol(XY)-1, ncol = length(lambdaSL))
                  for (i in 1:length(lambdaSL)) {
                    fit.UHuberLasso[,i] <- huberfit(XY[,-1],XY[,1], betam0 = beta.USqLasso, 
                                                    regularizer = "lasso", sig = sig, 
                                                    lambda = lambdaSL[i], itermax = 100)$betam
                  }
                  BIC <- apply(fit.UHuberLasso, 2, function(x) {
                    x[abs(x) < 1e-3] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberLasso <- fit.UHuberLasso[,which.min(BIC)]
                  beta.UHuberLasso[abs(beta.UHuberLasso) < 1e-3] <- 0
                  
                  ### MCP
                  fit.UHuberMCP <- matrix(NA, nrow = ncol(XY)-1, ncol = length(lambdaSL))
                  for (i in 1:length(lambdaSL)) {
                    fit.UHuberMCP[,i] <- huberfit(XY[,-1],XY[,1], betam0 = beta.USqLasso, 
                                                  regularizer = "MCP", sig = sig, 
                                                  lambda = lambdaSL[i], itermax = 100)$betam
                  }
                  BIC <- apply(fit.UHuberMCP, 2, function(x) {
                    x[abs(x) < 1e-3] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberMCP <- fit.UHuberMCP[,which.min(BIC)]
                  beta.UHuberMCP[abs(beta.UHuberMCP) < 1e-3] <- 0
                  
                  ### SCAD
                  fit.UHuberSCAD <- matrix(NA, nrow = ncol(XY)-1, ncol = length(lambdaSL))
                  for (i in 1:length(lambdaSL)) {
                    fit.UHuberSCAD[,i] <- huberfit(XY[,-1],XY[,1], betam0 = beta.USqLasso, 
                                                   regularizer = "SCAD", sig = sig, 
                                                   lambda = lambdaSL[i], itermax = 100)$betam
                  }
                  BIC <- apply(fit.UHuberSCAD, 2, function(x) {
                    x[abs(x) < 1e-3] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberSCAD <- fit.UHuberSCAD[,which.min(BIC)]
                  beta.UHuberSCAD[abs(beta.UHuberSCAD) < 1e-3] <- 0
                  
                  ### Lasso
                  return(list(USqLasso = round(beta.USqLasso,2),
                              USqMCP = round(beta.USqMCP,2),
                              USqSCAD = round(beta.USqSCAD,2),
                              UHuberLasso = round(beta.UHuberLasso,2),
                              UHuberMCP = round(beta.UHuberMCP,2),
                              UHuberSCAD = round(beta.UHuberSCAD,2)
                  ))
                })
  
  beta.Ufit <- NULL
  for (i in 1:6) {
    beta.Ufit <- cbind(beta.Ufit,c(sapply(fit,"[[",i)))  
  }
  eval.Ufit <- t(apply(beta.Ufit, 2, function(x) evallong(x, betam)))
  row.names(eval.Ufit) <- c("USqLasso",'USqMCP','USqSCAD',
                            "UHuberLasso",'UHuberMCP','UHuberSCAD')
  
  return(list(beta = beta.Ufit,
              eval = eval.Ufit
  ))
}

### fit only one huber lasso
Unitfit2 <- function(data) {
  fit <- lapply(data, 
                FUN = function(XY) {
                  np <- dim(XY)
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
                  fit.UHuberLasso <- hqreg(XY[,-(1:2)], XY[,1], lambda = lambdaSL)
                  BIC <- apply(fit.UHuberLasso$beta, 2, function(x) {
                    x[abs(x) < 1e-4] <- 0
                    log(sum((XY[,1]-XY[,-1] %*% x)^2)) + 
                      sum(x != 0) * log(np[1]) * log(np[2]) / (2 * np[1])
                  })
                  beta.UHuberLasso <- fit.UHuberLasso$beta[,which.min(BIC)]
                  beta.UHuberLasso[abs(beta.UHuberLasso) < 1e-4] <- 0
                  
                  
                  return(list(USqLasso = round(beta.USqLasso,2),
                              USqMCP = round(beta.USqMCP,2),
                              USqSCAD = round(beta.USqSCAD,2),
                              UHuberLasso = round(beta.UHuberLasso,2),
                              Huber_sig = fit.UHuberLasso$gamma
                  ))
                })
  
  beta.Ufit <- NULL
  for (i in 1:4) {
    beta.Ufit <- cbind(beta.Ufit,c(sapply(fit,"[[",i)))  
  }
  eval.Ufit <- t(apply(beta.Ufit, 2, function(x) evallong(x, betam)))
  row.names(eval.Ufit) <- c("USqLasso",'USqMCP','USqSCAD',
                            "UHuberLasso")
  
  return(list(beta = beta.Ufit,
              eval = eval.Ufit,
              sig = matrix(c(sapply(fit,"[[",5)), nrow = M, byrow = TRUE)
  ))
}



##### Fitting Huber Regression #####
huberfit <- function(XX, yy, regularizer = "lasso", omega0,
                     sig, lambda, itermax = 500) {
  
  lagrg <- huberloss <- object <- numeric()
  diff = iter = 1
  
  while (iter <= itermax & diff > 1e-3) {
    
    nu <- rep(0,length(betam0))
    eta <- threshold(regularizer, betam0 + nu/rho, lambda, rho, iota)
    nu <- nu + rho * (betam0 - eta)
    
    resid <- yy - XX %*% matrix(betam0, ncol=1)
    huberloss[iter] <- sum(ifelse(abs(resid) <= sig, 
                                  resid^2/2, sig*abs(resid)-sig^2/2))/np[1]
    object[iter] <- huberloss[iter] + 
      regularization(regularizer, betam0, lambda, iota)
    lagrg[iter] <- object[iter] + 
      sum(nu*(betam0 - eta)) + rho/2*sum((betam0 - eta)^2) 
    
    grad <- -t(ifelse(abs(resid) <= sig, resid, sign(resid)*sig)) %*% XX
    
    betam <- tau/(rho*tau+1) * (betam0/tau - grad/np[1] + rho*eta - nu)
    
    resid <- yy - XX %*% matrix(betam0, ncol=1)
    huberloss <- sum(ifelse(abs(resid) <= sig, 
                            resid^2/2, sig*abs(resid)-sig^2/2))/np[1]
    huberloss + sum(grad * (betam-betam0))/np[1] + sum((betam-betam0)^2)/2/tau +
      sum(nu*(betam - eta)) + rho/2*sum((betam - eta)^2) 
    huberloss + sum(nu*(betam0 - eta)) + rho/2*sum((betam0 - eta)^2) 
    
    betam[abs(betam) < 1e-4] <- 0
    eta[abs(eta) < 1e-4] <- 0
    
    diff <- norm(betam - betam0, "2")
    iter <- iter + 1
    betam0 <- betam
    rho <- rho * 1.03
  }
  return(list(object = object, lagrgian = lagrg, 
              betam = betam, iter = iter))
}

##### FedAvg #####
fedavg <- function(data, method = "Huber", regularizer = "lasso", omega0, sig, 
                   lambda, proportion, itermax = 200) {
  
  betaavg0 <- apply(matrix(omega0, ncol=M),1,mean)
  nu <- rep(0,length(betaavg0))
  lagrg <- huberloss <- object <- numeric()
  diff = iter = 1
  
  while (iter <= itermax & diff > 1e-3) {
    
    eta <- threshold(regularizer, betaavg0 + nu/rho, lambda, rho, iota)
    nu <- nu + rho * (betaavg0 - eta)
    
    client_sample <- sort(sample(M, M * proportion))
    grad <- rep(0, length(omega0))
    
    for (m in client_sample) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- betaavg0
      resid <- yy - XX %*% beta
      
      if(method == 'Huber') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(ifelse(abs(resid) <= sig[m], 
                                                 resid, sign(resid)*sig[m])) %*% XX
      }
      if(method == 'l2') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(resid) %*% XX
      }
    }
    
    gradavg <- apply(matrix(grad, ncol=M),1,sum)/proportion
    
    betaavg <- tau/(rho*tau+1) * (betaavg0/tau - gradavg/np[1]/M + rho*eta - nu)
    
    diff <- norm(betaavg - betaavg0, "2")
    iter <- iter + 1
    betaavg0 <- betaavg
  }
  return(list(beta = round(betaavg,2), iter = iter))
}