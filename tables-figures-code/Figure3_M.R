library(ggplot2)
library(dplyr)

Ms_list <- c(10,20,50,100)

# set path to the code folder
# path = setwd("..") 

##### plot for normal distribution #####

tab_M <- NULL

for (pp in 1:length(Ms_list)) {
  load(paste0('simulation-results/normal_M',Ms_list[pp],'_n100_p15.Rdata'))
  
  tab_M <- rbind(tab_M, cbind(apply(simplify2array(results), 1:2, median),
                              SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,1]),
                              lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1],
                              upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1])[,c(1,7,11:14)])
}
tab_M <- as.data.frame(cbind(M = rep(Ms_list,each=10),tab_M))
tab_M$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),length(Ms_list))

pdf(file = 'figures/normal_M_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_M %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')), 
       aes(x=M, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_line(linewidth=1.5)+geom_point(size=3) +
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method), 
              linetype=2, alpha=0.1) +
  xlab('Number of clients M') +
  ylab('log(MSE)') +
  ggtitle("Normal distribution") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  theme(legend.position=c(0.97,0), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()



##### plot for t distribution #####
tab_M <- NULL

for (pp in 1:length(Ms_list)) {
  load(paste0('simulation-results/t_M',Ms_list[pp],'_n100_p15.Rdata'))
  
  tab_M <- rbind(tab_M, cbind(apply(simplify2array(results), 1:2, median),
                              SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,1]),
                              lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1],
                              upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1])[,c(1,7,11:14)])
}
tab_M <- as.data.frame(cbind(M = rep(Ms_list,each=10),tab_M))
tab_M$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),length(Ms_list))

pdf(file = 'figures/t_M_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_M %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')), 
       aes(x=M, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_line(linewidth=1.5)+geom_point(size=3) +
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method), 
              linetype=2, alpha=0.1) +
  xlab('Number of clients M') +
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
  theme(legend.position=c(0.97,0), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()


##### plot for cauchy distribution #####

tab_M <- NULL

for (pp in 1:length(Ms_list)) {
  load(paste0('simulation-results/cauchy_M',Ms_list[pp],'_n100_p15.Rdata'))
  
  tab_M <- rbind(tab_M, cbind(apply(simplify2array(results), 1:2, median),
                              SD = sqrt(apply(simplify2array(results), 1:2, var, na.rm = TRUE)[,1]),
                              lower = apply(simplify2array(results), 1:2, quantile,0.25, na.rm = TRUE)[,1],
                              upper = apply(simplify2array(results), 1:2, quantile,0.75, na.rm = TRUE)[,1])[,c(1,7,11:14)])
}
tab_M <- as.data.frame(cbind(M = rep(Ms_list,each=10),tab_M))
tab_M$Method <- rep(c('Individual-Lasso','Individual-MCP','Individual-SCAD',
                      'Individual-Huber-Lasso',
                      'PerFL-l2-Lasso','PerFL-l2-MCP','PerFL-l2-SCAD',
                      'PerFL-RSR-Lasso','PerFL-RSR-MCP','PerFL-RSR-SCAD'),length(Ms_list))

pdf(file = 'figures/cauchy_M_mse.pdf', width = 9, height = 7.5)
ggplot(data = tab_M %>% 
         filter(Method %in% c('PerFL-RSR-MCP','PerFL-RSR-SCAD','PerFL-l2-MCP','PerFL-l2-SCAD')), 
       aes(x=M, y=log(MSE), colour=Method, linetype=Method)) + 
  geom_line(linewidth=1.5)+geom_point(size=3) +
  geom_ribbon(aes(ymin=log(lower), ymax=log(upper), fill = Method), 
              linetype=2, alpha=0.1) +
  xlab('Number of clients M') +
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
  theme(legend.position=c(0.97,0), legend.justification=c(1,0),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()



