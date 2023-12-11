##### Threshold Function #####
# thresholding operator for SCAD, MCP and LASSO regularizers
# td is the input vector
# lambda and iota are tuning parameters of the regularizer
# rho is the penalty parameter in the augmented Lagrangian

threshold <- function(regularizer = "MCP", td, lambda, rho, iota) {
  temp <- rep(NA, length(td))
  
  ### Soft threshold
  ST_td <- ifelse(abs(td) <= lambda/rho, 0, sign(td)*(abs(td)-lambda/rho))
  
  if(regularizer == 'MCP') {
    temp <- ifelse(abs(td) <= lambda*iota, ST_td/(1-1/iota/rho), td)
  } else if(regularizer == "SCAD") {
    temp[which(abs(td) <= lambda + lambda/rho)] <- ST_td[which(abs(td) <= lambda + lambda/rho)]
    temp[which(abs(td) > iota*lambda)] <- td[which(abs(td) > iota*lambda)]
    
    ST_td2 <- ifelse(abs(td) <= iota*lambda/rho/(iota-1), 0, 
                     sign(td)*(abs(td)-iota*lambda/rho/(iota-1))) / (1-1/(iota-1)/rho)
    temp[which((abs(td) <= iota*lambda) & (abs(td) > (lambda + lambda/rho)))] <- 
      ST_td2[which((abs(td) <= iota*lambda) & (abs(td) > (lambda + lambda/rho)))]
  } else if(regularizer == "lasso" | regularizer == "Lasso") {
    temp <- ST_td
  }
  return(temp)
}

###### Regularizers #####
# regularizers for SCAD, MCP and LASSO regularizers
# td is the input vector
# lambda and iota are tuning parameters of the regularizer

regularization <- function(regularizer = "MCP", td, lambda, iota) {
  temp <- rep(NA, length(td))
  if(regularizer == 'MCP') {
    temp <- ifelse(abs(td) <= lambda*iota, lambda*abs(td)-td^2/2/iota, lambda^2*iota/2)
  } else if(regularizer == "SCAD") {
    temp <- rep(NA, length(td))
    temp[which(abs(td) <= lambda)] <- (lambda*abs(td))[which(abs(td) <= lambda)]
    temp[which(abs(td) > iota*lambda)] <- lambda^2*(iota+1)/2
    
    td2 <- (iota*lambda*abs(td) - td^2/2-lambda^2/2)/(iota-1)
    temp[which((abs(td) > lambda) & (abs(td) <=  (iota*lambda)))] <- 
      td2[which((abs(td) > lambda) & (abs(td) <=  (iota*lambda)))]
  } else if(regularizer == "lasso" | regularizer == "Lasso") {
    temp <- lambda * abs(td)
  }
  return(sum(temp))
}