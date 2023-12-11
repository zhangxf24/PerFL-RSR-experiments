### This simulation takes about two hours for 
### one machine with eight cores

rm(list = ls(all = TRUE))

library(ncvreg)
library(ILAMM)
library(Matrix)
library(dplyr)
library(doParallel)

path = setwd("/Users/xfzhang/Documents/Research/HuberReg/Submissions/JASA submission/Revision/PerFL-RSR-experiments")

# set path to the code folder
# path = setwd("..") 

source('real-data-code/county_data_processing.R')
source('real-data-code/real_func.R')
source('simulation-code/Functions_regularizer.R')

##### setup #####
M <- length(county_list)

i <- rep(1:(M-1),(M-1):1)
j <- numeric()
for (m in 2:M) j <- c(j,m:M)
nM <- M*(M-1)/2
E <- sparseMatrix(i=c(1:nM,1:nM),j=c(i,j),
                  x=c(rep(1,nM),rep(-1,nM)))
OMEGA <- kronecker(E,diag(1,100))
ATA <- t(OMEGA) %*% OMEGA

rho <- 2
tau <- 1
iota <- 3 

r = (1+M)*rho*tau + 1
H <- (r-rho*tau)*diag(1, nrow = nrow(ATA)) - rho*tau*ATA

##### Model fitting function #####
crime_fun <- function(M) {
  ##### Assign training and testing data #####
  assigment <- train_test(crime_df)
  data <- assigment$data
  county_list <- assigment$county
  testDf <- as.data.frame(assigment$testDf)
  
  ##### Fitting by each unit #####
  Ufit <- Unitfit_real(data)
  evalUfit1 <- eval_real(round(Ufit$beta[,1],2), testDf)
  evalUfit2 <- eval_real(round(Ufit$beta[,2],2), testDf)
  
  ##### FedAvg #####
  lambdaSL = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))
  
  fit2_avg <- fedavg(data, method = "Huber", regularizer = "MCP", 
                     omega0 = Ufit$beta[,1], sig = rep(0.25,M), 
                     lambda = lambdaSL[6], proportion=1, itermax = 1000)
  fit3_avg <- fedavg(data, method = "Huber", regularizer = "SCAD",
                     omega0 = Ufit$beta[,1], sig = rep(0.25,M),
                     lambda = lambdaSL[6], proportion=1, itermax = 1000)
  
  evalavg2 <- eval_real(rep(fit2_avg$beta, M), testDf)
  evalavg3 <- eval_real(rep(fit3_avg$beta, M), testDf)
  
  ##### Proposed method #####
  
  fit2_pm <- PMfit(data, method = "Huber", regularizer = "MCP",
                   lambda1=lambdaSL[2], lambda2=lambdaSL[6], proportion = 1,
                   omega0 = Ufit$beta[,1], sig=rep(3,M), itermax = 1000)
  
  fit3_pm <- PMfit(data, method = "Huber", regularizer = "SCAD",
                   lambda1=lambdaSL[2], lambda2=lambdaSL[6], proportion = 1,
                   omega0 = Ufit$beta[,1], sig=rep(3,M), itermax = 1000)
  
  evalpm2 <- eval_real(as.vector(fit2_pm$omega), testDf)
  evalpm3 <- eval_real(as.vector(fit3_pm$omega), testDf)

  ##### evaluation #####
  
  eval_crime <-  list(evalUfit1, evalUfit2,
                      evalavg2, evalavg3, 
                      evalpm2, evalpm3, 
                      iter2 = fit2_pm$iter,
                      iter3 = fit3_pm$iter)
  
  return(eval_crime)
}

##### Simulation #####

m <- 10
set.seed(2023)

cl <- makeCluster(8)
registerDoParallel(cores=8)

results <- foreach(i = 1:m, .errorhandling = 'pass', 
               .packages = c('ncvreg','ILAMM')) %dopar% crime_fun(M)

stopCluster(cl)

save.image('simulation-results/crime_results_counties.Rdata')

