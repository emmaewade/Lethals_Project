#################
## Figure 2 and 4 -> Compare expected and inferred DFEs
#################

###################
####Prepare Panel 2
###################
plot_discrete_DFE <- function(alpha, beta, Nanc, type = "Inferred", prop = 0.0){
  # compute bin masses  
  x0 = (0. * 2 * Nanc)
  x1 = (1e-5 * 2 * Nanc) 
  x2 = (1e-4 * 2 * Nanc) 
  x3 = (1e-3 * 2 * Nanc) 
  x4 = (1e-2 * 2 * Nanc) 
  x5 = (1 * 2 * Nanc)
  
  print(x0)
  print(x1)
  print(x2)
  print(x3)
  print(x4)
  print(x5)
  
  if (prop == 0.1){
    print(prop)
  }
  
  a = pgamma(x1,shape=alpha,scale=beta) 
  a = a * (1 - prop)
  b = pgamma(x2,shape=alpha,scale=beta) - pgamma(x1,shape=alpha,scale=beta)
  b = b * (1 - prop)
  c = pgamma(x3,shape=alpha,scale=beta) - pgamma(x2,shape=alpha,scale=beta) 
  c = c * (1 - prop)
  d = pgamma(x4,shape=alpha,scale=beta) - pgamma(x3,shape=alpha,scale=beta) 
  d = d * (1 - prop)
  ef = 1 - pgamma(x4,shape=alpha,scale=beta) 
  e = ef * (1 - prop) + prop
  
  
  two_Nanc_s_range = c("Neutral", "Nearly Neutral", 
                       "Slightly Del.", "Moderately Del.", "Strongly Del.")
  
  
  check = tryCatch({
    
    df = data.frame("s" = two_Nanc_s_range, "Probability_mass" = c(a,b,c,d,e))
    
  }, error = function(e){
    
    print("Error!")
    
  })
  
  colnames(df) <- c("s", "Probability_mass")
  level_order <- two_Nanc_s_range

  #print(df)
  return(df)
}


# prepare data
com = read.table("Lethals_Project-main/Lethals_Project-main/figures/data/all_inferences.csv", h=T, sep = ",") 

#just need gamma distribution inferences
com$alpha = as.numeric(com$alpha)
com$beta = as.numeric(com$beta)
com$Nanc = as.numeric(com$Nanc)


ss = c("10", "100", "1000")
levels = c("0", "0.01", "0.05", "0.1")
h_co = c("0", "0.5")
rep = c(1:20)


df = data.frame()

#bins for inferred dfe

for (l in levels){
  
  for (h in h_co){
    
    for (s in ss) {
      
      for (r in rep){
        
        gam = com[which(com$dist == "gamma" & com$level == l & com$h_co == h & com$ss == s & com$rep == r), ]
        alpha_avg = gam$alpha
        beta_avg = gam$beta
        nanc_avg = gam$Nanc
        
        tmp = plot_discrete_DFE(alpha_avg, beta_avg, nanc_avg)
        tmp$type = "Inferred"
        tmp$level = l
        tmp$ss = s
        tmp$h = h
        tmp$rep = r
        
        check = tryCatch({
          
          df = rbind(df, tmp)
          
        }, error = function(e){
          
          print("Error!")
          
        })
        
      }
      
    }
  }
}


##bin expected dfe
for (l in levels ){
  for (h in h_co) {
    for (s in ss) {
      for (r in rep) {
        
        #take inferred 0%
        gam = com[which(com$dist == "gamma" & com$level == "0" & com$h_co == h & com$ss == s & com$rep == r), ]
        alpha_avg = gam$alpha
        beta_avg = gam$beta
        nanc_avg = gam$Nanc
        
        tmp_adj = plot_discrete_DFE(alpha_avg, beta_avg, nanc_avg, type = "Expected", prop = as.numeric(l))
        tmp_adj$type = "Expected"
        tmp_adj$level = l
        tmp_adj$ss = s
        tmp_adj$h = h
        tmp_adj$rep = r
        
        check = tryCatch({
          
          df = rbind(df, tmp_adj)
          
        }, error = function(e){
          
          print("Error!")
          
        })
        
      }
    }
  }
}

two_Nanc_s_range = c("Neutral", "Nearly Neutral", 
                     "Slightly Del.", "Moderately Del.", "Strongly Del.")
type = c("Expected", "Inferred")

##average replicates
to_plot = data.frame()
for (l in levels){
  for (h in h_co) {
    for (s in ss) {
      for (n in two_Nanc_s_range) {
        for (t in type) {
          tmp = df[which(df$ss == s & df$s == n & df$level == l & df$h == h & df$type == t), ]
          sd = sd(as.numeric(tmp$Probability_mass))
          mn = mean(as.numeric(tmp$Probability_mass))
          to_add = data.frame("ss" = s, "s"= n, "level" = paste(as.numeric(l)*100, "% Lethal", sep = ""), "h"= h, "sd" = sd, "mean" = mn, "Type" = t)
          to_plot = rbind(to_plot, to_add)
        }
      }      
    }
  } 
}

###############Prepare Panel 1 -- Dot Plot###########
#####################################################

#prepare data
df = read.table("Lethals_Project-main/Lethals_Project-main/figures/data/all_inferences.csv", h=T, sep = ",")

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

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

#h_co can be 0.5 or 0

library(ggplot2)
#to_plot$h can be 0 or 0.5
#to_plot$ss can be 10, 100, 1000

for (h in h_co){
  
  gg = ggplot(to_plot[which(to_plot$h == h & to_plot$ss == "1000"), ], aes(fill = factor (Type), x = factor(s, c("Neutral", "Nearly Neutral", 
                                                                                                                 "Slightly Del.", "Moderately Del.", "Strongly Del.")) , y = mean))+
    theme_bw()+
    geom_bar(position = "dodge", stat = "identity")+
    labs(fill = "DFE", x = "s", y = "Probability mass", fill = "Level")+
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20), plot.tag = element_text(), 
          legend.position = "right", axis.text.x=element_text(angle=25, vjust=1, hjust=1))+
    geom_errorbar(aes(ymin= mean  - sd, ymax= mean + sd), position=position_dodge(.9))+
    facet_wrap(~factor(level, levels=c("0% Lethal", "1% Lethal", "5% Lethal", "10% Lethal")), nrow = 4)+
    scale_fill_manual(values = c("grey75", "grey25"))
  

  gg2 = ggplot(df[which(df$dist == "gamma" & df$ss == "1000" & df$h_co == h), ], aes(alpha, beta, colour = factor(level)))+ 
    theme_bw()+
    geom_point(size = 5)+geom_point(aes(x=shape, y=kim_beta), colour="black", size = 5)+labs(colour="Lethal", x = "\u03B1", y = "\u03B2")+
    annotate("text", x = shape, y = kim_beta - 10, label = "\nKim 2017", color = "black", angle=0, size=8)+
    theme(text = element_text(size = 20), plot.tag = element_text(), legend.position = "right",  legend.background = element_rect(fill = "transparent"))+
    guides(color = guide_legend(keywidth = unit(0.0005, "inches"),
                                keyheight = unit(0.0005, "inches")))+
    scale_colour_manual(values = safe_colorblind_palette)
  

  library(cowplot)
  ggc = plot_grid(gg2, gg, nrow = 1, rel_widths = c(2/5, 3/5), labels = c('A', 'B'), label_size = 20)
  
  pdf(paste("figures/", "h_", h, "_exp_v_inferred_dfe.pdf", sep=""), width = 15, height = 8)
  print(ggc)
  dev.off()
  
}
