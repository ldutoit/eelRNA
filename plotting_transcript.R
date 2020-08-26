#This short script  plot and output read counts for transcripts.  

geneofinterest <- " TRINITY_DN2157_c0"  ## Enter atranscript or a prefix to get all the isoforms plotted (i.e. TRINITY_DN701_c2_g1_i2 or TRINITY_DN701_c2)

####MAIN 

library(ggplot2) # You can install ggplot2 by running install.packages("tidyverse")
## Define plottig function
plot_isoforms<-function(geneofinterest){
  sub_counts<-counts[grep("TRINITY_DN701_",rownames(counts)),]
  type = as.factor(rep(c("Silver","Yellow"),each=6))
  for (i in 1:dim(sub_counts)[1]){ # go through all the isoforms
    print(sub_counts[i,])
    transcript_name<-(rownames(sub_counts[i,]))
    if(transcript_name%in%rownames(results)){                     
      title_bit<-paste(transcript_name,", p-adjusted =",round(      results$adj.P.Val[which(rownames(results)==transcript_name)],4))}else{title_bit<-paste(transcript_name,"NOT passing filtering")}
    transcript<-data.frame(Counts= as.numeric(sub_counts[i,]),Type =type,color=color)
    x<-ggplot(data = transcript, mapping = aes(x = Type, y = Counts)) +  
      geom_jitter(alpha = 0.9, height=0,width=0.2) + theme_light() +ggtitle(title_bit)
    plot(x)
  }
}

load("countsand_logFCforplotting.RData") #Read counts and results are save in DE_analysis_els.Rmd
plot_isoforms(geneofinterest)