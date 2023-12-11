##### Fitting Evaluation #####
evalFit <- function(beta, beta0) {
  beta = matrix(beta, ncol = 1)
  beta0 = matrix(beta0, ncol = 1)
  s <- sum(beta0 != 0)
  p <- length(beta)
  beta = matrix(beta, ncol = 1)
  beta0 = matrix(beta0, ncol = 1)
  
  TP_RMean = sum(which(beta != 0) %in% which(beta0 != 0))
  FN_RMean = sum(which(beta == 0) %in% which(beta0 != 0))
  FP_RMean = sum(!which(beta != 0) %in% which(beta0 != 0))
  recall_RMean = TP_RMean / (TP_RMean + FN_RMean)
  precision_RMean = TP_RMean / (TP_RMean + FP_RMean)
  F1 <- 2*TP_RMean/(2*TP_RMean + FP_RMean + FN_RMean)
  CP <- all(which(beta0 != 0) %in% which(beta != 0))
  PCZ <- sum((beta == 0) * (beta0 == 0))/(p-s) 
  PICZ <- sum((beta == 0) * (beta0 != 0))/s
  res <-
    c(TP = TP_RMean/s,
      FN = FN_RMean/s,
      FP = FP_RMean/(length(beta)-s),
      recall = recall_RMean,
      precision = precision_RMean,
      F1 = F1,
      CP = CP,
      PCZ = PCZ,
      PICZ = PICZ)
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
             AMS = length(unique(round(omega,3))))
  eval
}