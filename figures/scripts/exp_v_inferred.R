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


### below is combining the 20 replicates, I assume I could probably condense this part
com = read.table("C:/big_summer/1_31/all_inference.csv", h=T, sep = ",") ### change path???
com_2 = read.table("C:/big_summer/1_31/second/all_inference.csv", h = T, sep = ",") ### change path??? 
com_2$rep = com_2$rep + 10
com = rbind(com, com_2)

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

library(ggplot)
#to_plot$h can be 0 or 0.5
#to_plot$ss can be 10, 100, 1000
gg = ggplot(to_plot[which(to_plot$h == "0.5" & to_plot$ss == "1000"), ], aes(fill = factor (Type), x = factor(s, c("Neutral", "Nearly Neutral", 
                                                       "Slightly Del.", "Moderately Del.", "Strongly Del.")) , y = mean))+
  theme_bw()+
  geom_bar(position = "dodge", stat = "identity")+
  labs(fill = "DFE", x = "s", y = "Probability mass", fill = "Level")+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 30), plot.tag = element_text(), 
        axis.text.x=element_text(size=15), legend.position = "right")+
  geom_errorbar(aes(ymin= mean  - sd, ymax= mean + sd), position=position_dodge(.9))+
  facet_wrap(~factor(level, levels=c("0% Lethal", "1% Lethal", "5% Lethal", "10% Lethal")), nrow = 4)

gg

