library(dplyr)
library(ggplot2)
library(tidyr)

# set path to the code folder
# path = setwd("..") 

load('simulation-results/crime_results.Rdata')

Client_MSE <- list()
for (i in 1:length(results)) {
  eval_temp <- results[[i]]
  temp <- NULL
  for (l in 1:6) {
    temp <- cbind(temp, eval_temp[[l]]$mse)
  }
  Client_MSE <- temp
}
Client_MSE <- apply(simplify2array(Client_MSE), 1:2, mean, na.rm=TRUE)

colnames(Client_MSE) <- c('Individual-MCP','Individual-SCAD',
                          'FedAvg-MCP','FedAvg-SCAD',
                          'PerFL-RSR-MCP','PerFL-RSR-SCAD')
tab_MSE <- Client_MSE %>% as.data.frame() %>%  
  gather(Method,MSE) %>% as.data.frame()

tab_MSE <- tab_MSE %>%
  mutate(Method = rep(c('Individual','FedAvg','PerFL-RSR'), each = M*2),
         Regularizer = rep(c('MCP','SCAD','MCP','SCAD','MCP','SCAD'), each = M))

pdf(file = 'figures/crime.pdf', width = 9, height = 7.5)
ggplot(data = tab_MSE %>% filter(MSE < 5), aes(x = Method, y = MSE, fill = Regularizer)) + 
  geom_boxplot(outlier.size=3,lwd=1) + 
  scale_color_brewer(palette = "Dark2") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",size=2,lineend = 10),
        axis.ticks.x=element_line(color="black",size=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1),
        legend.text=element_text(size=30),
        legend.title=element_text(size=30)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm")) 

dev.off()


