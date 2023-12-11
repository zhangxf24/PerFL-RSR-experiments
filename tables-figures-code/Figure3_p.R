library(ggplot2)
library(dplyr)

p <- c(15,100,200,300)

# set path to the code folder
# path = setwd("..") 

##### plot for normal distribution #####
tab_p <- NULL
for (t in 1:length(p)) {
  load(paste0('simulation-results/normal_M10_n100_p',p[t],'.Rdata'))
  
  tab_p <- rbind(tab_p, 
                 cbind(apply(simplify2array(results), 1:2, median),
                       SD = sqrt(apply(simplify2array(results), 1:2, var)[,1])/sqrt(50),
                       lower = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1],
                       upper = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1])[,c(1,7,11,12:14)])
}
tab_p <- as.data.frame(cbind(p = rep(p,each=10),tab_p))
tab_p$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),4)

pdf(file = 'figures/normal_p_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_p %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')),
       aes(x=p, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_point() + geom_line() + 
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method),
              linetype=2, alpha=0.1) +
  xlab('Covariate dimension p') +
  ylab('log(MSE)') +
  ggtitle("Normal distribution") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", linewidth=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  geom_line(size=1.5)+geom_point(size=3) +
  theme(legend.position=c(0.97,0), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()




##### plot for t distribution #####
tab_p <- NULL
for (t in 1:length(p)) {
  load(paste0('simulation-results/t_M10_n100_p',p[t],'.Rdata'))
  
  tab_p <- rbind(tab_p, 
                 cbind(apply(simplify2array(results), 1:2, median),
                       SD = sqrt(apply(simplify2array(results), 1:2, var)[,1]),
                       lower = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1],
                       upper = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1])[,c(1,7,11,12:14)])
}
tab_p <- as.data.frame(cbind(p = rep(p,each=10),tab_p))
tab_p$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),4)

pdf(file = 'figures/t_p_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_p %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')),
       aes(x=p, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_point() + geom_line() + 
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method), 
              linetype=2, alpha=0.1) +
  xlab('Covariate dimension p') +
  ylab('log(MSE)') +
  ggtitle("Student's t distribution") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  geom_line(size=1.5)+geom_point(size=3) +
  theme(legend.position=c(0.97,0), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()



##### plot for cauchy distribution #####
tab_p <- NULL
for (t in 1:length(p)) {
  load(paste0('simulation-results/cauchy_M10_n100_p',p[t],'.Rdata'))
  
  tab_p <- rbind(tab_p, 
                 cbind(apply(simplify2array(results), 1:2, median),
                       SD = sqrt(apply(simplify2array(results), 1:2, var)[,1]),
                       lower = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1],
                       upper = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1])[,c(1,7,11,12:14)])
}
tab_p <- as.data.frame(cbind(p = rep(p,each=10),tab_p))
tab_p$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),4)

pdf(file = 'figures/cauchy_p_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_p %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')),
       aes(x=p, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_point() + geom_line() + 
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method),
              linetype=2, alpha=0.1) +
  xlab('Covariate dimension p') +
  ylab('log(MSE)') +
  ggtitle("Cauchy distribution") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  geom_line(size=1.5)+geom_point(size=3) +
  theme(legend.position=c(0.52,0.62), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()



