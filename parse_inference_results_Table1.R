#here is some R code to parse the simulation/infernece results to make Table 1:
#This is the cleaned version of the code (updated April 2, 2023):


#read in data:
data<-read.csv("~/Dropbox/Kirk_stuff/Chris/all_inferences.aug.2022.csv")



###First, get differnce in log-likelihood for 0% datasets. This is how we define the empirical critical value for the LRT for the differnet mdoels.

ss<-1000 #which sample size
level<-0 #which amount of lethals
#get the data that we're going to analyze (we want h=0.5 & h=0, for the p_let=0 case):
use<-data[ data$ss==ss & data$level==level ,]

#now pull out the list of log-likelihoods for the relevant models:
gamma_ll<-use[use$dist=="gamma",3]
neu_gamma_ll<-use[use$dist=="neugamma",3]
let_gamma_ll<-use[use$dist=="gammalet",3]
neu_let_gamma_ll<-use[use$dist=="neugammalet",3]

#now get LRT critical value for the lethal_gamma vs gamma model:
LRT_let<-2*(let_gamma_ll-gamma_ll)
LRT_let_use<-ifelse(LRT_let<0,0,LRT_let)
quantile(LRT_let_use,0.95)
#18.12781 #round to 20 to be conservative


#now get LRT critical value for the neutral_gamma vs gamma model:
LRT_neu<-2*(neu_gamma_ll-gamma_ll)
LRT_neu_use<-ifelse(LRT_neu<0,0,LRT_neu)
quantile(LRT_neu_use,0.95)
#4.391985 #round to 5 to be conservative


#now get LRT critical value for the neutral_gamma_lethal vs gamma model:
LRT_neu_let_gamma<-2*(neu_let_gamma_ll-gamma_ll)
LRT_neu_let_gamma_use<-ifelse(LRT_neu_let_gamma<0,0,LRT_neu_let_gamma)
quantile(LRT_neu_let_gamma_use,0.95)
#17.32003 #round to 20 to be conservative





#here is a function to determine whether the more complex models fit better than simpler models, when using the desired critical value. This is used to assess the number of datasets where more complex models fit better

calc_llike_diff<-function(i,ss,level,h)
{
	d_like<-c(0,0,0)
	use<-data[data$rep==i & data$ss==ss & data$level==level & data$h_co==h,]

	#now get log-likelihood of gamma
	gamma_ll<-use[use$dist=="gamma",3]

	#now get log likelihood of neut_gamma
	neu_gamma_ll<-use[use$dist=="neugamma",3]
	#now get log likelihood of lethal gamma
	let_gamma_ll<-use[use$dist=="gammalet",3]
	#now get log-likelihood of lethal+gamma+neutral
	neu_let_gamma_ll<-use[use$dist=="neugammalet",3]

	if(length(gamma_ll)>0)
		{
			d_like[1]<-2*(let_gamma_ll-gamma_ll)
			d_like[2]<-2*(neu_gamma_ll-gamma_ll)
			d_like[3]<-2*(neu_let_gamma_ll-gamma_ll)
		}
		
		return(d_like)
}

#now set up some loops to find the proportion of simulated datasets where the more compelx models fit better:

ss<-1000
out<-numeric()
#llike_thresh<-5 #use for Gamma vs. Neut_gamma, comment out for the others
llike_thresh<-20 #use for Gamma vs. Let_gamma as well as Gamma vs. Neut_Let_Gamma; comment out for the Gamma vs. Neut_gamma



for (h in c(0,0.5))
{	

	for (level in c(0,0.01,0.05,0.1))
	{
		
		#now loop over all 20 reps
		count_lethal<-0
		count_neu<-0
		count_neu_let<-0
		for (i in seq(1,20,by=1))
		{
			
			like_diff_out<-calc_llike_diff(i,ss,level,h)
			if(like_diff_out[1]>llike_thresh)
			{
				count_lethal<-count_lethal+1
			}
			
			if(like_diff_out[2]>llike_thresh)
			{
				count_neu<-count_neu+1
			}
			
			if(like_diff_out[3]>llike_thresh)
			{
				count_neu_let<-count_neu_let+1
			}
			
		}
		
		
		use2<-data[ data$ss==ss & data$level==level & data$h_co==h & data$dist=="gammalet", ]
		
		mean_lethal<-mean(as.numeric(use2$plet))
		
		print(c(h,ss,level,count_lethal,count_neu,mean_lethal))
		out_vec<-c(h,ss,level,mean_lethal,count_lethal,count_neu,count_neu_let,llike_thresh)
		out<-rbind(out,out_vec)
	
	}
	
	
}	
colnames(out)<-c("h","samp_size","true_p_let","est_p_let","gamma+let","gamma+neu","gamma+neu+let","llike_cutoff")

file_name <- paste("~/Dropbox/Kirk_stuff/Chris/inference_results.table.aug.2022.llike_cutoff_" , llike_thresh, ".txt")
write.table(out,file_name,quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)



 





