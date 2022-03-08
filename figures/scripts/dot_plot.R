#prepare data
df = read.table("all_inferences.csv", h=T, sep = ",")

#gamma inferences
df_gamma = df[which(df$dist == "gamma"), ]
df_gamma$alpha = as.numeric(df_gamma$alpha)
df_gamma$beta = as.numeric(df_gamma$beta)

#rescale
df$mean = (df$alpha * df$beta * 2) / (2 * df$Nanc)

df = df_gamma

#rescale Kim 2017 values
mean = 0.01314833
shape = 0.186
Nanc = 10000
scale = mean/shape

#beta = mean*(2*Nanc) / alpha * 2 = beta
kim_beta = mean*(2*Nanc) / (shape * 2)


require(utf8)
library(ggplot2)

gg2 = ggplot(df[which(df$dist == "gamma" & df$ss == "1000" & df$h_co == "0.5"), ], aes(alpha, beta, colour = factor(level)))+ 
  theme_bw()+
  geom_point(size = 5)+geom_point(aes(x=shape, y=kim_beta), colour="black", size = 5)+labs(colour="Lethal", x = "\u03B1", y = "\u03B2")+
  annotate("text", x = shape, y = kim_beta - 10, label = "\nKim 2018", color = "black", angle=0, size=8)+
  theme(text = element_text(size = 30), plot.tag = element_text(), legend.position = "right",  legend.background = element_rect(fill = "transparent"))+
  guides(color = guide_legend(keywidth = unit(0.0005, "inches"),
                              keyheight = unit(0.0005, "inches")))
gg2

