library(xtable)

##### Table 1 part 1 #####

# set path to the code folder
# path = setwd("..") 

load('simulation-results/BIC_M10_n50.Rdata')

eval <- list()
for(mm in 1:50) {
  eval[[mm]] <- results[[mm]]$eval
}

rows <- c(5:8,11:12)
cols <- c(1,5,6,7,11)

means <- apply(simplify2array(eval), 1:2, mean, na.rm=TRUE)[rows,cols]
sds <- apply(simplify2array(eval), 1:2, sd, na.rm=TRUE)[rows,cols]
means <- round(means,2)
sds <- round(sds/sqrt(50),2)
means
sds
xtable(cbind(paste0(means[,1],'(',sds[,1],')'),means[,2:4],
             paste0(means[,5],'(',sds[,5],')')))

##### Table 1 part 2 #####

load('simulation-results/BIC_M10_n100.Rdata')

eval <- list()
for(mm in 1:50) {
  eval[[mm]] <- results[[mm]]$eval
}

rows <- c(5:8,11:12)
cols <- c(1,5,6,7,11)

means <- apply(simplify2array(eval), 1:2, mean, na.rm=TRUE)[rows,cols]
sds <- apply(simplify2array(eval), 1:2, sd, na.rm=TRUE)[rows,cols]
means <- round(means,2)
sds <- round(sds/sqrt(50),2)
means
sds
xtable(cbind(paste0(means[,1],'(',sds[,1],')'),means[,2:4],
             paste0(means[,5],'(',sds[,5],')')))

##### Table 1 part 3 #####

load('simulation-results/BIC_M50_n100.Rdata')

eval <- list()
for(mm in 1:50) {
  eval[[mm]] <- results[[mm]]$eval
}

rows <- c(5:8,11:12)
cols <- c(1,5,6,7,11)

means <- apply(simplify2array(eval), 1:2, mean, na.rm=TRUE)[rows,cols]
sds <- apply(simplify2array(eval), 1:2, sd, na.rm=TRUE)[rows,cols]
means <- round(means,2)
sds <- round(sds/sqrt(50),2)
means
sds
xtable(cbind(paste0(means[,1],'(',sds[,1],')'),means[,2:4],
             paste0(means[,5],'(',sds[,5],')')))








