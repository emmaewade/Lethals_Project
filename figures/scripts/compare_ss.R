plot_discrete <- function(alpha, beta, Nanc, prop = 0){
  # compute bin masses  
  x0 = 0. * 2 * Nanc
  x1 = 1e-5 * 2 * Nanc
  x2 = 1e-4 * 2 * Nanc
  x3 = 1e-3 * 2 * Nanc
  x4 = 1e-2 * 2 * Nanc
  x5 = 1 * 2 * Nanc
  
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
    colnames(df) <- c("s", "Probability_mass")
    level_order <- two_Nanc_s_range
    print(df)
    return(df)
    
  }, error = function(e){
    
    return()
  
  })
  
}

#prepare data
com = read.table("c:/big_summer/1_31/all_inference.csv", h=T, sep=',', stringsAsFactors = F)
com_2 = read.table("c:/big_summer/1_31/second/all_inference.csv", h=T, sep=',', stringsAsFactors = F)
com_2$rep = com_2$rep + 10
com = rbind(com, com_2)
com$alpha = as.numeric(com$alpha)
com$beta = as.numeric(com$beta)
com$Nanc = as.numeric(com$Nanc)

ss = c("10", "100", "1000")
levels = c("0", "0.01", "0.05", "0.1")
h_co = c("0", "0.5")
rep = c(1:20)

df = data.frame()

##inferred 
for (l in levels){
  
  for (h in h_co){
    
    for (s in ss) {
      
      for (r in rep) {
      
        gam = com[which(com$dist == "gamma" & com$level == l & com$h_co == h & com$ss == s & com$rep == r), ]
        tmp = plot_discrete(gam$alpha, gam$beta, gam$Nanc)
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


###expected dfe
for (l in levels ){
  for (h in h_co) {
      for (r in rep) {
        
        #take 0% lethal dfe
        gam = com[which(com$dist == "gamma" & com$level == "0" & com$h_co == h & com$ss == "1000" & com$rep == r), ]
        tmp = plot_discrete(gam$alpha, gam$beta, gam$Nanc, as.numeric(l))
        tmp$level = l
        tmp$ss = "Expected"
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

two_Nanc_s_range = c("Neutral", "Nearly Neutral", 
                     "Slightly Del.", "Moderately Del.", "Strongly Del.")

types = c("10", "100", "1000", "Expected")

to_plot = data.frame()
for (l in levels){
  for (h in h_co) {
    for (s in types) {
      for (n in two_Nanc_s_range) {
        tmp = df[which(df$ss == s & df$s == n & df$level == l & df$h == h), ]
        sd = sd(as.numeric(tmp$Probability_mass))
        mn = mean(as.numeric(tmp$Probability_mass))
        to_add = data.frame("ss" = s, "s"= n,  "level" = paste(as.numeric(l)*100, "% Lethal", sep = ""), "h"= h, "sd" = sd, "mean" = mn)
        to_plot = rbind(to_plot, to_add)
      }      
    }
  } 
}



gg = ggplot(to_plot[which(to_plot$h == "0"), ], aes(fill = ss, x = factor(s,levels= c("Neutral", "Nearly Neutral", 
                                                                                                 "Slightly Del.", "Moderately Del.", "Strongly Del.")) , y = mean))+
  theme_bw()+
  geom_bar(position = "dodge", stat = "identity")+
  facet_wrap(~factor(level), nrow = 4)+
  labs( x = "s", y = "Mean Probability Mass", fill = "Sample Size")+
  geom_errorbar(aes(ymin= mean  - sd, ymax= mean + sd), position=position_dodge(.9))+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20), plot.tag = element_text(), 
        axis.text.x=element_text(size=13))


#####dot plot######################
df = data.frame()
ss = c("10", "100", "1000")
levels = c("0", "0.01", "0.05", "0.1")
h_co = c("0", "0.5")
rep = c(1:20)

for (l in levels){
  
  for (h in h_co){
    
    to_add = data.frame("ss" = "Expected", "level" = l, "h" = h, "sd" = 0, "mean" = as.numeric(l))
    df = rbind(df, to_add)
    
    for (s in ss) {
      
        gamma_let = com[which(com$dist == "gammalet" & com$level == l & com$h_co == h & com$ss == s), ]
        mean = mean(as.numeric(gamma_let$plet))
        sd = sd(as.numeric(gamma_let$plet))
        to_add = data.frame("ss" = s, "level" = l, "h"= h, "sd" = sd, "mean" = mean)
        
        check = tryCatch({
          
          df = rbind(df, to_add)
          
        }, error = function(e){
          
          print("Error!")
          
        })
      
    }
  }
}

#h can be 0 or 0.5
gg2 = ggplot(df[which(df$h == "0.5"), ], aes(x = factor(level), y = as.numeric(mean), colour = factor(ss, levels = c("10", "100", "1000", "Expected"))))+ 
  theme_bw()+
  #geom_errorbar(aes(ymin= mean  - sd, ymax= mean + sd))+
  geom_point(size = 5)+
  #facet_wrap(~factor(level), nrow = 4)+
  labs(x = "True Lethal Proportion", y = expression(paste("Mean ", hat("p")["let"])), colour="Sample Size")+
  theme(text = element_text(size = 20),
        plot.tag = element_text(), legend.position = "right",  
        legend.background = element_rect(fill = "transparent"))
  
gg2


library(cowplot)
ggc = plot_grid(gg, gg2, nrow = 1, rel_widths = c(3/5, 2/5), labels = c('A', 'B'), label_size = 20)
ggc
