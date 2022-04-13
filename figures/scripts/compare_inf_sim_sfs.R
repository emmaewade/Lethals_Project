
library(RColorBrewer)                            # Load RColorBrewer
library(ggplot2)

#returns sfs dataframe to plot
stacked_sfs <- function(del, sample, num_sites, titlename = NA) {
  
  #remove unneeded columns
  del_sfs = del[3:(sample+2)]
  del_sfs = t(del_sfs)
  del_tbl = data.frame("Sites" = 1:sample, "Number" = del_sfs, "Type" = "Deleterious")
  colnames(del_tbl) = c("Sites", "Number", "Type")
  
  tmp_df = rbind(del_tbl)
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

#################################
## prepare simulated SFS 
##########
sfs_10 = read.table("sim_avg_sfs10.csv", sep = ',', h = T)
sfs_100 = read.table("sim_avg_sfs100.csv", sep = ',', h = T)
sfs_1000 = read.table("sim_avg_sfs1000.csv", sep = ',', h = T)

prop = c("0", "0.01", "0.05", "0.1")
print_prop = c("Simulated - 0.00", "Simulated - 0.01", "Simulated - 0.05", "Simulated - 0.10")
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
    
    stacked_10 = stacked_sfs(df_10[which(df_10$type == "nonsyn"), ], 10, 10)
    stacked_10$ss = "10"
    stacked_10$h_co = h
    stacked_10$Source = print_prop[i]
    
    stacked_100 = stacked_sfs(df_100[which(df_100$type == "nonsyn"), ], 100, 10)
    stacked_100$ss = "100"
    stacked_100$h_co = h
    stacked_100$Source = print_prop[i]
    
    stacked_1000 = stacked_sfs(df_1000[which(df_1000$type == "nonsyn"), ], 1000, 10)
    stacked_1000$ss = "1000"
    stacked_1000$h_co = h
    stacked_1000$Source = print_prop[i]
    
    df = rbind(df, stacked_10, stacked_100, stacked_1000)
    
    
  }
  
}

sim_df = df

#################################
## prepare inferred SFS 
##########

sfs_10 = read.table("computed_avg_sfs10.csv", sep = ',', h = T)
sfs_100 = read.table("computed_avg_sfs100.csv", sep = ',', h = T)
sfs_1000 = read.table("computed_avg_sfs1000.csv", sep = ',', h = T)

prop = c("0", "0.01", "0.05", "0.1")
print_prop = c("Inferred - 0.00", "Inferred - 0.01", "Inferred - 0.05", "Inferred - 0.10")
h_co = c("0", "0.5")
dist = c("gamma")

i = 0

num_sites = 10
df = data.frame()

for (p in prop){
  
  i = i + 1
  
  for (h in h_co) {
    
    for (d in dist) {
    
      df_10 = sfs_10[which(sfs_10$prop == p & sfs_10$h_co == h & sfs_10$dist == d), ]
      df_100 = sfs_100[which(sfs_100$prop == p & sfs_100$h_co == h & sfs_100$dist == d), ]
      df_1000 = sfs_1000[which(sfs_1000$prop == p & sfs_1000$h_co == h & sfs_1000$dist == d), ]
      
      stacked_10 = stacked_sfs(df_10, 10, 10)
      stacked_10$ss = "10"
      stacked_10$h_co = h
      stacked_10$Source = print_prop[i]
      
      stacked_100 = stacked_sfs(df_100, 100, 10)
      stacked_100$ss = "100"
      stacked_100$h_co = h
      stacked_100$Source = print_prop[i]
      
      stacked_1000 = stacked_sfs(df_1000, 1000, 10)
      stacked_1000$ss = "1000"
      stacked_1000$h_co = h
      stacked_1000$Source = print_prop[i]
      
      df = rbind(df, stacked_10, stacked_100, stacked_1000)
    
    }
    
  }
  
}

#combine fitdadi inferred sfs and simulated sfs
df = rbind(df, sim_df)

#ss can be 10, 100, or 1000
gg = ggplot(df[which(df$ss == 1000),], aes(fill=factor(Source, levels = c("Inferred - 0.00", "Simulated - 0.00", "Inferred - 0.01", "Simulated - 0.01", "Inferred - 0.05", "Simulated - 0.05", "Inferred - 0.10", 'Simulated - 0.10')),
                    y=Number, x=Sites)) + 
  geom_bar(position='dodge', stat='identity')+
  theme_bw()+
  scale_x_discrete(limits = c(1:num_sites,"+"), position  = "bottom")+
  scale_color_brewer(palette = "Set2")+
  labs(fill = "")+
  ylab("Number of SNPs")+
  xlab("SNP Frequency")+
  theme(text = element_text(size = 30), legend.position = "bottom")+
  facet_wrap(~h_co)

gg
