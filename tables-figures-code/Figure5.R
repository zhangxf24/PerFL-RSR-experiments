library(ggplot2)
library(dplyr)
library(tidyr)

# set path to the code folder
# path = setwd("..") 

##### Figure for proportion #####
load('simulation-results/proportion.Rdata')

res <- list()
for (i in 1:m) {
  tab <- NULL
  for (j in 1:length(props)) {
    tab <- rbind(tab, cbind(proportion = props[j], results[[i]][[j]][,1:200]))
  }
  res[[i]] <- as.matrix(tab)
}

tab_mean <- apply(simplify2array(res), 1:2, mean)
tab_lower <- apply(simplify2array(res), 1:2, quantile, 0.25, na.rm = TRUE)
tab_upper <- apply(simplify2array(res), 1:2, quantile, 0.75, na.rm = TRUE)
tab_mean <- tab_mean %>% as.data.frame()
tab_lower <- tab_lower %>% as.data.frame()
tab_upper <- tab_upper %>% as.data.frame()

tab_mean$Method <- tab_lower$Method <- tab_upper$Method <- 
  rep(c('SqLasso','SqMCP','SqSCAD', 'PMLasso','PMMCP','PMSCAD'), length(props))
tab_mean$type <- rep('MSE', 6*length(props))
tab_lower$type <- rep('lower', 6*length(props))
tab_upper$type <- rep('upper', 6*length(props))

colnames(tab_mean) <- c('Proportion', 1:200, 'Method', 'type')
colnames(tab_lower) <- c('Proportion', 1:200, 'Method', 'type')
colnames(tab_upper) <- c('Proportion', 1:200, 'Method', 'type')

tab_mean <- tab_mean %>% 
  gather(Iteration,MSE,-c(Proportion,Method,type)) %>%
  filter(!is.na(MSE)) %>% 
  mutate(Iteration = as.numeric(Iteration),
         Proportion = as.factor(Proportion))
# tab_mean %>% glimpse

tab_lower <- tab_lower %>% 
  gather(Iteration,lower,-c(Proportion,Method,type)) %>%
  filter(!is.na(lower)) %>% 
  mutate(Iteration = as.numeric(Iteration),
         Proportion = as.factor(Proportion))
# tab_lower %>% glimpse

tab_upper <- tab_upper %>% 
  gather(Iteration,upper,-c(Proportion,Method,type)) %>%
  filter(!is.na(upper)) %>% 
  mutate(Iteration = as.numeric(Iteration),
         Proportion = as.factor(Proportion))
# tab_upper %>% glimpse

tab_prop <- tab_mean %>% mutate(communication = Iteration*M*as.numeric(Proportion))
tab_prop$lower <- tab_lower$lower
tab_prop$upper <- tab_upper$upper
# tab_prop %>% glimpse

pdf(file = 'figures/Proportion_vs_iter_MCP.pdf', width = 9, height = 7.5)
ggplot(data = tab_prop %>% 
         filter(Method=='PMMCP' & Proportion %in% (c(1,3,5,7,9)/10)) %>%
         filter(Iteration <= 200), 
       aes(x=Iteration, y=log(MSE), colour=Proportion, linetype=Proportion)) + 
  geom_line(size=1.5) +
  xlab('Iteration number T') +
  ylab('log(MSE)') +
  ggtitle('Regularizer: MCP') + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill=NA,color="black", size=3, linetype="solid"),
        axis.ticks.y=element_line(color="black",linewidth=2,lineend = 10),
        axis.ticks.x=element_line(color="black",linewidth=2,lineend = 10),
        axis.title.x = element_text(vjust=0.1)) +
  theme(text = element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1),
        legend.text=element_text(size=30),
        legend.title=element_text(size=30)) +
  theme(legend.background = element_rect(fill=rgb(1,1,1,alpha = 0.001), colour = NA)) +
  theme(plot.margin = margin(1,1,1,1, "cm"))

dev.off()
