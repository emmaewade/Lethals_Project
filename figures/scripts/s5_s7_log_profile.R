#read inference
create_tbl <- function(data_dir){
  
  ss = c("10", "100", "1000")
  rep = as.factor(1:20)
  h_co = c("0.0", "0.5")
  levels = c("0.00", "0.01", "0.05", "0.10")
  prop = c(as.character("0.0"), as.character(seq(0.02, 0.5, by= 0.02)))
  
  columns = c("level", "ss", "h_co", "log", "sim_plet",
              "plet", "pneu",
              "alpha", "beta", "rep")
  
  df = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(df) = columns
  
  for (s in ss) {
    
    for (l in levels){
      
      for (h in h_co) {
        
        for (p in prop){
          
          for (r in rep) {
            
            path = paste(data_dir, "DFE_inference_", p, "_", l, "_", h, "_", s, "_", r, ".csv", sep = "")
             
            tb = read.table(path, header = T, sep = ',') 
            
            tb$rep = r
            
            tb$sim_plet = l
            tb$h_co = h
            tb$ss = s
            
            df = rbind(tb, df)  
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(df)
}

#create plots
plot <- function(df, h) {
  
  library(ggplot2)
  
  
  props = c("0", "0.01", "0.05", "0.1")
  reps = c(1:20)
  
  
  graphs = list()
  
  for (p in props){
    
    t = data.frame()
    
    #subtract each replicate likelihood from highest likelihood
    for (r in reps) {
    
      temp = df[which(df$level == p & df$rep == r & df$h_co == h & df$ss == 1000),] 
      temp$Likelihood = as.numeric(temp$Likelihood) + as.numeric(-max(temp$Likelihood))
      
      t = rbind(temp, t)
    
    }
    
    #find max mle
    max = t$plet[which(t$Likelihood == 0)]
    max = as.numeric(max)
    mean_max = mean(max, na.rm=T)
    
    std_var = sd(max, na.rm=T)
    
    gg <- ggplot(t, aes(x=plet, y=Likelihood)) + 
      theme_bw()+
      geom_line(aes(col=factor(rep, levels = as.character(c(seq(1,20, by=1)))))) +
      labs(title = paste(as.numeric(p)*100, "% Lethals", sep = ""), colour = "Replicate: ", x = expression(hat("p")["let"]), y = "Log likelihood") + 
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5), text = element_text(size = 20))+
      geom_vline(xintercept=as.numeric(p), linetype='dashed', color = 'red')+
      geom_vline(xintercept=as.numeric(mean_max), linetype='dashed', color = 'blue')+
      annotate("text", x = as.numeric(p)+0.003, y = -20, label = expression(paste("Expected ", hat("p")["let"])), color = "red", angle=90, size=5) +
      annotate("text", x = as.numeric(mean_max) + 0.003, y = -10, label = "Mean MLE", color = "blue", angle=90, size=5) +
      annotate("text", x = .112, y = -2, label = paste("sd = ", round(std_var,5)), size=5) +
      coord_cartesian(xlim = c(-0.002, 0.13), ylim=c(-32,0), expand = F)+
      scale_colour_grey()
      
    graphs[[p]] = gg
    
  }
    
  return(graphs)
}

#prepare data
t_df = read.table("Lethals_Project-main/Lethals_Project-main/figures/data/all_log_20_inference.csv", h=T, sep=",")

#h 0 or 0.5
h_co = c(0, 0.5)
for (h in h_co) {
  gg = plot(t_df, h)
  
  require(gridExtra)
  require("ggpubr")
  ggc <- ggarrange(gg[["0"]], gg[["0.01"]], gg[["0.05"]], gg[["0.1"]], common.legend = FALSE, ncol=2, nrow=2)
  
  pdf(paste("figures/h_", h, "_log_profile.pdf", sep = ""), width = 15, height = 10)
  print(ggc)
  dev.off()
}

