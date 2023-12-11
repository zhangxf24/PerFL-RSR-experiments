###### Evaluation #####

eval_real <- function(betav, testDf) {
  size <- sum(apply(matrix(betav, ncol = M), 1, function(x) sum(unique(x) != 0))) + 1
  sse <- mse <- numeric()
  np <- length(betav)/M
  for (i in 1:length(county_list)) {
    temp <- testDf %>% filter(countystate == county_list[i])
    coef <- betav[((i-1)*np+1):(i*np)]
    # size[i] <- sum(coef != 0)
    XX <- as.matrix(temp[,c(101,1:99)])
    YY <- temp[,100]
    sse[i] <- sum((YY - XX %*% coef)^2)
    mse[i] <- mean((YY - XX %*% coef)^2)
  }
  return(list(size = size, sse = sse, mse = mse, county = county_list))
}

##### Fit by each state #####
Unitfit_real <- function(data) {
  fit <- lapply(data, 
                FUN = function(XY) {
                  np <- dim(XY)
                  XY <- as.matrix(XY)
                  lambdaSL <- max(abs(t(XY[,1]) %*% XY[,-1])) / np[1] *
                    exp(seq(log(10^1), log(10^-6), length = 20))
                  
                  ### Fitting by each unit with Huber loss ###
                  
                  ### Lasso
                  fit.UHuberMCP <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "MCP", 
                                                  lSeq = lambdaSL, nfolds = 2)
                  fit.UHuberSCAD <- cvNcvxHuberReg(XY[,-(1:2)], XY[,1], penalty = "SCAD", 
                                                   lSeq = lambdaSL, nfolds = 2)
                  
                  return(list(UHuberMCP = fit.UHuberMCP$beta,
                              UHuberSCAD = fit.UHuberSCAD$beta,
                              tau.min = c(fit.UHuberMCP$tauMin, fit.UHuberSCAD$tauMin)))
                })
  
  beta.Ufit <- NULL
  for (i in 1:2) {
    beta.Ufit <- cbind(beta.Ufit,c(sapply(fit,"[[",i)))  
  }
  
  
  
  return(list(beta = beta.Ufit,
         tau.min = matrix(c(sapply(fit,"[[",3)), nrow = length(data), byrow = TRUE)))
}

##### Proposed method #####

PMfit <- function(data, method = "Huber", regularizer = "lasso", omega0, sig,
                  lambda1, lambda2, proportion, itermax = 200) {
  
  lagrg <- huberloss <- object <- numeric()
  diff = iter = 1
  gamma <- rep(1,nrow(OMEGA))
  nu <- rep(1,length(omega0))
  
  while (iter <= itermax & diff > 1e-5) {
    
    delta <- threshold(regularizer, OMEGA %*% omega0 + gamma/rho, lambda1, rho, iota)
    eta <- threshold(regularizer, omega0 + nu/rho, lambda2, rho, iota)
    gamma <- gamma + rho * (OMEGA %*% omega0- delta)
    nu <- nu + rho * (omega0 - eta)
    
    huloss <- 0
    for (m in 1:M) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      np <- dim(XX)
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      huloss <- huloss +
        sum(ifelse(abs(resid) <= sig[m], resid^2/2, sig[m]*abs(resid)-sig[m]^2/2))/length(yy)
    }
    huberloss[iter] <- huloss/M
    object[iter] <- huberloss[iter] + regularization(regularizer, delta, lambda1, iota) +
      regularization(regularizer, eta, lambda2, iota)
    lagrg[iter] <- object[iter] +
      sum(gamma*(OMEGA %*% omega0 - delta)) + rho/2*sum((OMEGA %*% omega0 - delta)^2) +
      sum(nu*(omega0 - eta)) + rho/2*sum((omega0 - eta)^2)
    
    Delta_omega <- H %*% omega0 / r -
      tau / r * ( - rho*t(OMEGA) %*% (delta - gamma/rho) - rho*eta + nu)
    
    client_sample <- sort(sample(M, M * proportion))
    grad <- rep(0, length(omega0))
    for (m in client_sample) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      np <- dim(XX)
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      
      if(method == 'Huber') {
        grad[(1:np[2])+(m-1)*np[2]] <- (-t(ifelse(abs(resid) <= sig[m],
                                                  resid, sign(resid)*sig[m])) %*% XX) / length(yy)
      }
      if(method == 'l2') {
        grad[(1:np[2])+(m-1)*np[2]] <- (-t(resid) %*% XX)/length(yy)
      }
      
    }
    
    omega <- Delta_omega - tau / r * grad/M
    
    diff <- norm(omega - omega0, "2")
    iter <- iter + 1
    omega0 <- omega
  }
  omega <- round(omega,2)
  return(list(object = object, lagrgian = lagrg, omega = omega, iter = iter))
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
      np <- dim(XX)
      beta <- betaavg0
      resid <- yy - XX %*% beta

      if(method == 'Huber') {
        grad[(1:np[2])+(m-1)*np[2]] <- (-t(ifelse(abs(resid) <= sig[m],
                                                  resid, sign(resid)*sig[m])) %*% XX)/length(yy)
      }
      if(method == 'l2') {
        grad[(1:np[2])+(m-1)*np[2]] <- (-t(resid) %*% XX)/length(yy)
      }
    }

    gradavg <- apply(matrix(grad, ncol=M),1,sum)/proportion

    betaavg <- tau/(rho*tau+1) * (betaavg0/tau - gradavg/M + rho*eta - nu)

    diff <- norm(betaavg - betaavg0, "2")
    iter <- iter + 1
    betaavg0 <- betaavg
  }
  return(list(beta = round(betaavg,2), iter = iter))
}

##### Assign training and testing data #####

train_test <- function(crime_df) {

  data <- list()
  testDf <- NULL
  
  for (i in 1:length(county_list)) {
    temp <- crime_df %>% filter(countystate == county_list[i])
    sam <- sample(nrow(temp), 0.8*nrow(temp))
    
    data[[i]] <- as.matrix(temp[sam,c(100,101,1:99)])
    testDf <- rbind(testDf, temp[-sam,]) 
  }
  
  return(list(data = data, testDf = testDf, county = county_list))
}




