
##### Proposed Method #####

PMfit <- function(data, method = "Huber", regularizer = "lasso", omega0, sig, 
                  lambda1, lambda2, proportion, itermax = 200) {
  
  # initialization
  lagrg <- huberloss <- object <- numeric()
  diff = iter = 1
  gamma <- rep(0,nrow(OMEGA))
  nu <- rep(0,length(omega0))
  
  while (iter <= itermax & diff > 1e-3) {
    
    # update parameters other than model weights
    delta <- threshold(regularizer, OMEGA %*% omega0 + gamma/rho, lambda1, rho, iota)
    eta <- threshold(regularizer, omega0 + nu/rho, lambda2, rho, iota)
    gamma <- gamma + rho * (OMEGA %*% omega0- delta)
    nu <- nu + rho * (omega0 - eta)
    
    # calculate the vaules of objective function and Lagrangian function
    huloss <- 0
    for (m in 1:M) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      huloss <- huloss+sum(ifelse(abs(resid) <= sig[m], resid^2/2, sig[m]*abs(resid)-sig[m]^2/2))
    }
    huberloss[iter] <- huloss/M/np[1]
    object[iter] <- huberloss[iter] + regularization(regularizer, delta, lambda1, iota) +
      regularization(regularizer, eta, lambda2, iota)
    lagrg[iter] <- object[iter] + 
      sum(gamma*(OMEGA %*% omega0 - delta)) + rho/2*sum((OMEGA %*% omega0 - delta)^2) +
      sum(nu*(omega0 - eta)) + rho/2*sum((omega0 - eta)^2) 
    
    # calcultate the transfer information from server to clients
    Delta_omega <- H %*% omega0 / r -
      tau / r * ( - rho*t(OMEGA) %*% (delta - gamma/rho) - rho*eta + nu)
    
    # sample clients
    client_sample <- sort(sample(M, M * proportion))
    
    # calculate gradient on each client
    grad <- rep(0, length(omega0))
    for (m in client_sample) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      
      if(method == 'Huber') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(ifelse(abs(resid) <= sig[m], 
                                                 resid, sign(resid)*sig[m])) %*% XX
      }
      if(method == 'l2') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(resid) %*% XX
      }
      
    }
    
    # update model weights
    omega <- Delta_omega - tau / r * grad/M/np[1] 
    
    diff <- norm(omega - omega0, "2")
    iter <- iter + 1
    omega0 <- omega
  }
  omega <- round(omega,2)
  return(list(object = object, lagrgian = lagrg, omega = omega, iter = iter))
}

PMfit_simu <- function(data, method = "Huber", regularizer = "lasso", omega0, sig, 
                       lambda1, lambda2, proportion, itermax = 200) {
  
  mse <- lagrg <- huberloss <- object <- numeric()
  diff = iter = 1
  
  gamma <- rep(0,nrow(OMEGA))
  nu <- rep(0,length(omega0))
  
  while (iter <= itermax & diff > 1e-3) {
    
    mse[iter] <- norm(omega0-c(betam), type = "2")
    
    delta <- threshold(regularizer, OMEGA %*% omega0 + gamma/rho, lambda1, rho, iota)
    eta <- threshold(regularizer, omega0 + nu/rho, lambda2, rho, iota)
    gamma <- gamma + rho * (OMEGA %*% omega0- delta)
    nu <- nu + rho * (omega0 - eta)
    
    huloss <- 0
    for (m in 1:M) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      huloss <- huloss+sum(ifelse(abs(resid) <= sig[m], resid^2/2, sig[m]*abs(resid)-sig[m]^2/2))
    }
    huberloss[iter] <- huloss/M/np[1]
    object[iter] <- huberloss[iter] + regularization(regularizer, delta, lambda1, iota) +
      regularization(regularizer, eta, lambda2, iota)
    lagrg[iter] <- object[iter] + 
      sum(gamma*(OMEGA %*% omega0 - delta)) + rho/2*sum((OMEGA %*% omega0 - delta)^2) +
      sum(nu*(omega0 - eta)) + rho/2*sum((omega0 - eta)^2) 
    
    Delta_omega <- H %*% omega0 / r -
      tau / r * ( - rho*t(OMEGA) %*% (delta - gamma/rho) - rho*eta + nu)
    
    ### sample clients
    client_sample <- sort(sample(M, M * proportion))
    grad <- rep(0, length(omega0))
    for (m in client_sample) {
      yy <- data[[m]][,1]
      XX <- data[[m]][,-1]
      beta <- omega0[(1:np[2])+(m-1)*np[2]]
      resid <- yy - XX %*% beta
      
      if(method == 'Huber') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(ifelse(abs(resid) <= sig[m], 
                                                 resid, sign(resid)*sig[m])) %*% XX
      }
      if(method == 'l2') {
        grad[(1:np[2])+(m-1)*np[2]] <- -t(resid) %*% XX
      }
      
    }
    
    omega <- Delta_omega - tau / r * grad/M/np[1] 
    
    diff <- norm(omega - omega0, "2")
    iter <- iter + 1
    omega0 <- omega
  }
  omega <- round(omega,2)
  return(list(mse = mse, object = object, lagrgian = lagrg, omega = omega, iter = iter))
}