library(dplyr)
library(xtable)

# set path to the code folder
# path = setwd("..") 

load('simulation-results/crime_results_counties.Rdata')

##### Table 2 #####

n1 <- NULL
for (i in 1:length(county_list)) {
  temp <- crime_df %>% filter(countystate == county_list[i])
  n1 <- c(n1, nrow(temp))
}
summary(n1)

##### Table 3 #####

MSEs <- sapply(results, function(res) sapply(res[1:6], function(x) mean(x$mse)))
MSEs_sd <- sapply(results, function(res) sapply(res[1:6], function(x) sd(x$mse)))
Sizes <- sapply(results, function(res) sapply(res[1:6], function(x) x$size))

means <- rbind(apply(t(MSEs), 2, mean),
             apply(t(MSEs_sd), 2, mean),
             apply(t(Sizes), 2, mean))
sds <- rbind(apply(t(MSEs), 2, sd)/sqrt(10),
             apply(t(MSEs_sd), 2, mean)/sqrt(10),
             apply(t(Sizes), 2, sd)/sqrt(10))
means <- t(round(means, 2))
sds <- t(round(sds, 2))
means
sds
xtable(matrix(paste0(means[,],'(',sds[,],')'), nrow = 6))

