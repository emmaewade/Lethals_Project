##########
## stacked sfs
##########

library(ggplot2)

#returns sfs dataframe to plot
stacked_sfs <- function(del, neu, leth, sample, num_sites, titlename = NA) {
  
  #remove unneeded columns
  neu_sfs = neu[3:(sample+2)]
  neu_sfs = t(neu_sfs)
  neu_tbl = data.frame("Sites" = 1:sample, "Number" = neu_sfs, "Type" = "Neutral")
  colnames(neu_tbl) = c("Sites", "Number", "Type")
  
  #remove unneeded columns
  let_sfs = leth[3:(sample+2)]
  let_sfs = t(let_sfs)
  let_tbl = data.frame("Sites" = 1:sample, "Number" = let_sfs, "Type" = "Lethal")
  colnames(let_tbl) = c("Sites", "Number", "Type")
  
  #remove unneeded columns
  del_sfs = del[3:(sample+2)]
  del_sfs = t(del_sfs)
  del_sfs = del_sfs - let_sfs
  del_tbl = data.frame("Sites" = 1:sample, "Number" = del_sfs, "Type" = "Deleterious")
  colnames(del_tbl) = c("Sites", "Number", "Type")
  
  tmp_df = rbind(del_tbl, neu_tbl, let_tbl)
  colnames(tmp_df) = c("Sites", "Number", "Type")
  
  temp_df = data.frame()
  #colnames(df) = c("Sites", "Number", "Type")
  t = data.frame()
  
  #find frequency between num_sites and length
  for (ty in unique(tmp_df$Type)){
    
    sumhigh = sum(tmp_df$Number[which(tmp_df$Sites > num_sites & tmp_df$Type == ty)])
    
    cond_df = data.frame("Sites" = 1:num_sites, "Number" = tmp_df$Number[which((tmp_df$Sites <= num_sites) & tmp_df$Type == ty)], "Type" = ty)
    
    high_df = data.frame("Sites" = "+", "Number" = sumhigh, "Type" = ty)
    
    t = rbind(cond_df, high_df)
    temp_df = rbind(temp_df, cond_df, high_df)
    
  }
  
  
  return(temp_df)
  
}

sfs_10 = read.table("sim_avg_sfs10.csv", sep = ',', h = T)
sfs_100 = read.table("sim_avg_sfs100.csv", sep = ',', h = T)
sfs_1000 = read.table("sim_avg_sfs1000.csv", sep = ',', h = T)

prop = c("0", "0.01", "0.05", "0.1")
print_prop = c("0% Lethal", "1% Lethal", "5% Lethal", "10% Lethal")
h_co = c("0", "0.5")

i = 0

num_sites = 10
df = data.frame()

for (p in prop){
  
  i = i + 1
  
  for (h in h_co) {
    
    df_10 = sfs_10[which(sfs_10$prop == p & sfs_10$h_co == h), ]
    df_100 = sfs_100[which(sfs_100$prop == p & sfs_100$h_co == h), ]
    df_1000 = sfs_1000[which(sfs_1000$prop == p & sfs_1000$h_co == h), ]
    
    stacked_10 = stacked_sfs(df_10[which(df_10$type == "nonsyn"), ], df_10[which(df_10$type == "syn"), ], df_10[which(df_10$type == "lethal"), ], 10, 10)
    stacked_10$ss = "10"
    stacked_10$h_co = h
    stacked_10$prop = print_prop[i]
    
    stacked_100 = stacked_sfs(df_100[which(df_100$type == "nonsyn"), ], df_100[which(df_100$type == "syn"), ], df_100[which(df_100$type == "lethal"), ], 100, 10)
    stacked_100$ss = "100"
    stacked_100$h_co = h
    stacked_100$prop = print_prop[i]
    
    stacked_1000 = stacked_sfs(df_1000[which(df_1000$type == "nonsyn"), ], df_1000[which(df_1000$type == "syn"), ], df_1000[which(df_1000$type == "lethal"), ], 1000, 10)
    stacked_1000$ss = "1000"
    stacked_1000$h_co = h
    stacked_1000$prop = print_prop[i]
    
    df = rbind(df, stacked_10, stacked_100, stacked_1000)
    
    
  }
  
}

df$Sites[which(df$Sites == "+")] = num_sites + 1
barwidth = 0.35
num_sites = 10

sample = 1000 #can be 10, 100, or 1000
hco = 0.5 #can be 0 or 0.5

ggplot() + 
  theme_bw()+
  geom_bar(data = df[which(df$Type == "Neutral" & df$ss == sample & df$h_co == hco), ],
           mapping = aes(x = Sites, y = as.numeric(Number), fill = Type), 
           position='stack', stat='identity', width = barwidth)+
  geom_bar(data = df[which(df$ss == sample & df$h_co == hco & (df$Type == "Deleterious" | df$Type == "Lethal")), ], 
           mapping = aes(x = as.numeric(Sites) + barwidth + 0.01, y = as.numeric(Number), fill = Type), 
           position='stack', stat='identity', width = barwidth)+
  scale_x_discrete("SNP Frequency", limits = c(1:(num_sites+1)), labels = c(1:num_sites, "+"), position  = "bottom")+
  labs(y = "Number of SNPs")+
  theme(plot.title = element_text(hjust = 0.5),  text = element_text(size = 12))+
  facet_wrap(~factor(prop, levels = c("0% Lethal", "1% Lethal", "5% Lethal", "10% Lethal")))

gg

########################################
#### 1% Number of Lethal
lethal_1 = df[which(df$prop == "1% Lethal" & df$Type == "Lethal"), ]

#facet = simulated level
gg = ggplot(lethal_1, aes(fill=ss,
                    y=Number, x=Sites)) + 
  geom_bar(position='dodge', stat='identity')+
  theme_bw()+
  scale_x_discrete(limits = c(1:num_sites,"+"), position  = "bottom")+
  scale_color_brewer(palette = "Set2")+
  labs(fill = "Sample Size")+
  ylab("Number of lethals")+
  xlab("Lethal frequency")+
  theme(text = element_text(size = 13), legend.position = "right")+
  facet_wrap(~h_co)

gg
