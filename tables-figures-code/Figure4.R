library(ggplot2)
library(dplyr)

# set path to the code folder
# path = setwd("..") 

##### Figure for Huber sigma #####

sigmal <- 1:7

tab_sigma <- NULL

for (t in 1:length(sigmal)) {
  load(paste0('simulation-results/lambda2_',t,'.Rdata'))
  
  tab_sigma <- rbind(tab_sigma, cbind(apply(simplify2array(results), 1:2, median),
                                      SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,1]),
                                      lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1],
                                      upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1])[,c(1,7,11:14)])
}

tab_sigma <- as.data.frame(cbind(sigma = rep(sigmal,each=3),tab_sigma))
tab_sigma$Method <- rep(c('PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'), length(sigmal))

pdf(file = 'figures/sigma_MSE.pdf', width = 9, height = 7.5)
ggplot(data = tab_sigma, 
       aes(x=sigma, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_point() + geom_line() +
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method), 
              linetype=2, alpha=0.1) +
  # xlab(expression("Privacy budget " ~ epsilon)) +
  xlab(expression('Huber loss parameter' ~ sigma)) +
  ylab('log(MSE)') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  # ggtitle("FapLAD")+
  theme(text = element_text(size=30)) +
  geom_line(linewidth=1.5)+geom_point(size=3) +
  theme(legend.position=c(0.97,1), legend.justification=c(1,1),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()




##### Figure for lambda1#####

lambdaSL1 = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))

tab_lb <- NULL

for (t in 1:(length(lambdaSL1))) {
  load(paste0('simulation-results/lambda1_',t,'.Rdata'))
  
  tab_lb <- rbind(tab_lb, cbind(apply(simplify2array(results), 1:2, mean),
                                  SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,11]),
                                  lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,11],
                                  upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,11])[,c(1,7,11:14)])
}
tab_lb <- as.data.frame(cbind(lambda = rep(lambdaSL1,each=3),tab_lb))
tab_lb$Method <- rep(c('PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'), length(lambdaSL1))

pdf(file = 'figures/lambda1_AMS.pdf', width = 9, height = 7.5)
ggplot(data = tab_lb, 
       aes(x=log(lambda), y=AMS, colour=Method, linetype=Method)) + 
  geom_line(linewidth=1.5)+geom_point(size=3) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Method), 
              linetype=2, alpha=0.1) +
  xlab(expression('Regularization parameter log(' ~ lambda[1] ~')')) +
  ylab('AMS') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  theme(legend.position=c(0.98,1), legend.justification=c(1,1),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()


##### Figure for lambda2#####

lambdaSL1 = exp(seq(log(10 ^ 1), log(10 ^ -6), length = 10))

tab_lb <- NULL

for (t in 1:(length(lambdaSL1))) {
  load(paste0('simulation-results/lambda2_',t,'.Rdata'))
  
  tab_lb <- rbind(tab_lb, cbind(apply(simplify2array(results), 1:2, mean),
                                SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,7]),
                                lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,7],
                                upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,7])[,c(1,7,11:14)])
}
tab_lb <- as.data.frame(cbind(lambda = rep(lambdaSL1,each=3),tab_lb))
tab_lb$Method <- rep(c('PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'), length(lambdaSL1))

pdf(file = 'figures/lambda2_F1.pdf', width = 9, height = 7.5)
ggplot(data = tab_lb, 
       aes(x=log(lambda), y=F1, colour=Method, linetype=Method)) + 
  geom_line(linewidth=1.5)+geom_point(size=3) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = Method), 
              linetype=2, alpha=0.1) +
  xlab(expression('Regularization parameter log(' ~ lambda[2] ~')')) +
  ylab(expression( ~F[1]~ 'Score')) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  theme(legend.position=c(0.52,1), legend.justification=c(1,1),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()




