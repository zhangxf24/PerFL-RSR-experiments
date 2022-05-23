### Fitting Evaluation
evalFit <- function(beta, beta0) {
  beta = matrix(beta, ncol = 1)
  beta0 = matrix(beta0, ncol = 1)
  s <- sum(beta0 != 0)
  p <- length(beta)
  beta = matrix(beta, ncol = 1)
  beta0 = matrix(beta0, ncol = 1)
  # MSE = sum((beta - beta0)^2)
  TP_RMean = sum(which(beta != 0) %in% which(beta0 != 0))/s
  FN_RMean = sum(which(beta == 0) %in% which(beta0 != 0))/s
  FP_RMean = sum(!which(beta != 0) %in% which(beta0 != 0))/s
  recall_RMean = TP_RMean / (TP_RMean + FN_RMean)
  precision_RMean = TP_RMean / (TP_RMean + FP_RMean)
  F1 <- 2/(1/recall_RMean+1/precision_RMean)
  CP <- all(which(beta0 != 0) %in% which(beta != 0))
  PCZ <- sum((beta == 0) * (beta0 == 0))/(p-s) 
  PICZ <- sum((beta == 0) * (beta0 != 0))/s
  res <-
    c(TP = TP_RMean,
      FN = FN_RMean,
      FP = FP_RMean,
      recall = recall_RMean,
      precision = precision_RMean,
      F1 = F1,
      CP = CP,
      PCZ = PCZ,
      PICZ = PICZ)
  # colnames(res) <- c('l2', 'TP', 'FN', 'FP', 'recall', 'precision')
  res
}

evallong <- function(omega, betam) {
  M <- ncol(betam)
  eval <- rep(0,9)
  for (m in 1:M) {
    eval <- eval + evalFit(omega[(1:np[2])+(m-1)*np[2]], betam[,m])
  }
  eval  <- c(MSE = norm(omega-c(betam), type = "2"), 
                      eval/M,
                      AMS = length(unique(omega)))
  eval
}
  

### Fitting by each unit with square loss

UnitSq <- function(data, penalty, betam) {
  fit <- lapply(data,
                FUN = function(XY) {
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  fit <- ncvreg(XY[,-(1:2)], XY[,1],
                                penalty = penalty, lambda = lambdaSL)
                  
                  BIC <- apply(fit$beta, 2, function(x) {
                    np[1]*log(sum((XY[,1]-XY[,-1] %*% x)^2)/np[1]) + 
                      sum(x != 0) * log(np[1]) 
                  })
                  beta <- fit$beta[,which.min(BIC)]
                  beta[abs(beta) < 10^-4] = 0
                  beta
                })
  eval <- rep(0,9)
  M <- length(data)
  for (m in 1:M) {
    eval <- eval + evalFit(fit[[m]], betam[,m])
  }
  eval <- c(MSE = norm(unlist(fit)-c(betam), type = "2"), 
            eval/M,
            AMS = length(unique(unlist(fit))))
  return(list(eval=eval, omega = unlist(fit)))
}

### Fitting by each unit with Huber loss

UnitHu <- function(data, penalty, betam, tau = 10) {
  fit <- lapply(data,
                FUN = function(XY) {
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  BIC <- 100000
                  beta <- rep(0,np[2])
                  for (i in 1:length(lambdaSL)) {
                    fit0 <- ncvxHuberReg(XY[,-1], XY[,1], lambda = lambdaSL[i], 
                                         penalty = penalty, tau = tau,
                                         iteMax = 500, itcpIncluded = TRUE) 
                    betat <- fit0$beta
                    betat[abs(betat) < 10^-4] = 0
                    resid <- XY[,1]-XY[,-1] %*% betat
                    huber.loss <- sum(ifelse(abs(resid) <= tau, resid^2/2, tau*abs(resid)-tau^2/2))
                    BICt <- np[1] * log(huber.loss/np[1]) +  sum(betat != 0) * log(np[1])
                    if(BICt < BIC) {
                      beta <- betat
                      BIC <- BICt
                    }
                  }
                  beta
                })
  eval <- rep(0,9)
  M <- length(data)
  for (m in 1:M) {
    eval <- eval + evalFit(fit[[m]], betam[,m])
  }
  eval <- c(MSE = norm(unlist(fit)-c(betam), type = "2"), 
            eval/M,
            AMS = length(unique(unlist(fit))))
  return(eval)
}

UnitHu.tf <- function(data, penalty, betam) {
  fit <- lapply(data,
                FUN = function(XY) {
                  fit0 <- tfNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = penalty)
                  fit0$beta
                })
  eval <- rep(0,9)
  M <- length(data)
  for (m in 1:M) {
    eval <- eval + evalFit(fit[[m]], betam[,m])
  }
  eval <- c(MSE = norm(unlist(fit)-c(betam), type = "2"), 
            eval/M,
            AMS = length(unique(unlist(fit))))
  return(eval)
}


### proposed

PHHS_Huber <- function(data, penalty = "MCP", betam, omega, tau = 100, lambda = 1, eta=1, itermax = 100) {
  M <- length(data)
  i <- rep(1:(M-1),(M-1):1)
  j <- numeric()
  for (m in 2:M) j <- c(j,m:M)
  nM <- M*(M-1)/2
  E <- sparseMatrix(i=c(1:nM,1:nM),j=c(i,j),
                    x=c(rep(1,nM),rep(-1,nM)))
  A <- rbind(kronecker(E,diag(1,np[2])), diag(1,M*np[2]))
  ATA <- t(A) %*% A
  phi.ATA <- eigen(ATA)$values
  rho <- (1+sqrt(33))/min(phi.ATA)
  
  nu <- ifelse(penalty == "MCP", 3, 3.7)
  
  lambda <- 1
  eta <- 1
 
  r = max(phi.ATA)*rho*eta+1/M
  H <- r*diag(1, nrow = M*np[2]) - rho*eta*ATA
  # range(eigen(H)$values)
  
  if(missing(omega)) {
    omega <- rnorm(M*np[2],1)
    delta <- A %*% omega
    gamma <- rep(0,nrow(A))
  }
  
  iter = 0
  mse <- 1
  while (iter <= itermax & mse > 0.0001) {
    td <- A %*% omega + gamma/rho
    if(penalty == 'MCP') {
      delta <- ifelse(abs(td) <= lambda*nu, 
                      ifelse(td <= lambda/rho, 0, sign(td)*(abs(td)-lambda/rho))/(1-1/nu/rho),td)
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
      delta <- ifelse(td <= lambda/rho, 0, sign(td)*(abs(td)-lambda/rho))
    }
    
    grad <- numeric()
    for (m in 1:M) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- omega[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      grad[(1:np[2])+(m-1)*np[2]] <- -t(ifelse(abs(resid) <= tau, resid, sign(resid)*tau)) %*% XX
    }
    
    omega <- H %*% omega / r - 
      eta / r * (grad/M/np[1] - t(A) %*% (rho * delta - gamma))
    
    gamma <- gamma + rho*(A %*% omega - delta)
    
    iter <- iter + 1
    mse <- norm(omega-c(betam), type = "2")
  }
  eval <- rep(0,9)
  M <- length(data)
  omega <- as.vector(omega)
  for (m in 1:M) {
    eval <- eval + evalFit(omega[(1:np[2])+(m-1)*np[2]], betam[,m])
  }
  eval <- c(MSE = norm(omega-c(betam), type = "2"), 
            eval/M,
            AMS = length(unique(omega)))
  return(eval)
}











