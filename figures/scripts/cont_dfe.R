#prepare data
com = read.table("c:/big_summer/repo/all_inferences.csv", h=T, sep=',', stringsAsFactors = F)
com$alpha = as.numeric(com$alpha)
com$beta = as.numeric(com$beta)
com$Nanc = as.numeric(com$Nanc)

ss = c("10", "100", "1000")
levels = c("0", "0.01", "0.05", "0.1")
h_co = c("0", "0.5")


df = data.frame()

#for each inference
for (l in levels){
  
  for (h in h_co){
    
    for (s in ss) {
      
      gam = com[which(com$dist == "gamma" & com$level == l & com$h_co == h & com$ss == s), ]
      alpha = mean(gam$alpha)
      beta = mean(gam$beta)
      nanc = mean(gam$Nanc)
      mean= (alpha * beta * 2) / (2 * nanc)#0.01266228
      
      #rescale
      scale=mean/alpha
      
      g <- rgamma(n=100000,shape = alpha, scale= scale)
      tmp <- data.frame("dist" = -g, "level" = l, "h" = paste("h =",h), "ss" = s)
      df = rbind(df, tmp)
      
    }
  }
}

library(ggplot2)

#ss - 10, 100, or 1000
ggplot(df[which(df$ss == "1000"),], aes(x = (dist + 1), fill = factor(level, levels = c("0.1", "0.05", "0.01", "0"))))+
  theme_bw()+
  geom_histogram(bins = 10000,position = "identity")+
  coord_cartesian(xlim = c(0.50, 1), ylim = c(0, 75))+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20))+
  labs(x = "s", y = "Frequency", fill = "Lethals")+
  facet_wrap(~factor(h))
