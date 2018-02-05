library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)



################Helper Functions#######################
rank_matrix_values<-function(matrix,margin=2){
  #Rank the expression level for a gene in each sample. Each column is a sample.
  ranked_matrix=apply(matrix,MARGIN=margin,FUN=rank)
  return(ranked_matrix)
}

report_quantiles<-function(v_data,probs=seq(0,1,0.1),type=7,na.rm=FALSE){
  #Input is a named vector of measurements
  #Output is a named vector of the quantiles for the measurements
  quantiles=unique(quantile(x=v_data,type=type,probs=probs,na.rm=na.rm))
  
  number_of_groups=length(quantiles)-1
  
  
  quantile_labled_data=cut(v_data,breaks=quantiles,include_lowest=TRUE,labels=1:number_of_groups)
  
  names(quantile_labled_data)=names(v_data)
  
  return(quantile_labled_data)
}

get_rows_coorsponding_to_samples<-function(df,unique_sample_string,ignore.case=TRUE,column_index=1){
  #Returns the row index of if unique_sample_string is found in the column specified by column_index
  #for that row
  cn=df[,column_index]
  return(grep(unique_sample_string,cn,ignore.case = ignore.case))
}



##The Function that performs the filtering specifically for the NF54 and PB58 Datasets
identify_dephasing_drivers_nf54_pb58<-function(m_fpkm_data,m_ref_data,
                                               sample1_name,sample2_name,
                                               minimum_acceptable_correlation,
                                               quantiles=seq(0,1,0.01),
                                               timepoint,
                                               rank_diff_outfile="rank_diff.csv",
                                               rank_diff_hist_outfile="rank_diff_histogram.pdf",
                                     abs_quantiles_outfile="abs_values_quantiles.csv",outdir="Correlation_Improvement_Graphs",
                                     randomize_data=FALSE,
                                     randomization_seed=5712){
  
#m_fpkm_data is a matrix of the experimental data
#m_ref_data is a matrix of the reference data used to check the life-cycle correlations
#sample1_name and sample2_name are the column names in m_fpkm_data for the samples whose correlations
  #are being computed
#minimum_acceptable_correlation is stopping criteria for the re-phasing algorithm
#quantiles is a vector indicating which quantiles to make
#timepoint is a string that can be used to filter the samples based on the timepoint indicated in the column name
  #and is very specific to this specific dataset
  
#If randomize_data is TRUE will randomize the ranks of the experimental da

  #Create the out directory
  dir.create(outdir)
  ##Get the ranks of the values

  #Get the absolute rank difference for each gene between the samples
  #If randomize_data is TRUE the abs_rank_diff for a gene will be randomized. Otherwise, the true rank difference will
  #be used
  if(isTRUE(randomize_data)){
    print("Randomizing Data")
    
    #Randomize the experimental data measurements
    set.seed(randomization_seed) #Setting seed to make randomization repeatable
    nrows=nrow(m_fpkm_data)
    ncols=ncol(m_fpkm_data)
    number_of_measurements=nrows*ncols
    matrix_mean=mean(m_fpkm_data)
    matrix_std=sd(m_fpkm_data)
    matrix_colnames=colnames(m_fpkm_data)
    matrix_rownames=rownames(m_fpkm_data)
    m_fpkm_data=matrix(rnorm(n=number_of_measurements,mean=matrix_mean,sd=matrix_std),nrows,ncols)
   
    
    colnames(m_fpkm_data)=matrix_colnames
    rownames(m_fpkm_data)=matrix_rownames
    View(m_fpkm_data)
    
    m_fpkm_ranked=rank_matrix_values(m_fpkm_data)
    v_sample1=m_fpkm_ranked[,sample1_name]
    v_sample2=m_fpkm_ranked[,sample2_name]
    
    
    
    rank_diff=v_sample2-v_sample1
    df_rankdiff=data.frame(gene=names(rank_diff),rank.diff=rank_diff)
    write.csv(df_rankdiff,rank_diff_outfile)
    abs_rank_diff=abs(rank_diff)
    
    #Plot out the rank differences as a histogram
    pdf(rank_diff_hist_outfile)
    hist(rank_diff,main="Randomized Rank Difference Distribution", xlab="Rank Difference")
    dev.off()
    
  }else{
    print("Not Randomizing Data")
    
    m_fpkm_ranked=rank_matrix_values(m_fpkm_data)
    v_sample1=m_fpkm_ranked[,sample1_name]
    v_sample2=m_fpkm_ranked[,sample2_name]
    
    rank_diff=v_sample2-v_sample1
    df_rankdiff=data.frame(gene=names(rank_diff),rank.diff=rank_diff)
    write.csv(df_rankdiff,rank_diff_outfile)
    abs_rank_diff=abs(rank_diff)
    
    #Plot out the rank differences as a histogram
    pdf(rank_diff_hist_outfile)
    hist(rank_diff,main="True Rank Difference Distribution", xlab="Rank Difference")
    dev.off()
  }
  
  

  #Get the quantiles for the rank differences
  abs_rank_quantiles=report_quantiles(abs_rank_diff,probs=quantiles,type=7,na.rm=TRUE)
  df_quantiles=data.frame(gene=names(abs_rank_quantiles),quantile=abs_rank_quantiles)
  write.csv(df_quantiles,abs_quantiles_outfile)
  quantiles=sort(unique(abs_rank_quantiles),decreasing=TRUE)
  
  old_ref_cor=cor(x=v_sample1,y=v_sample2,method="spearman")
  
  #Loop to perform the filtering
  i=1
  number_of_genes=nrow(df_fpkm)
  quantile_size=ceiling(number_of_genes/length(quantiles))
  #Continue while data correlation is less than minimum correlation and there is still unfiltered data
  while(old_ref_cor<minimum_acceptable_correlation & nrow(m_fpkm_data) >quantile_size){ 
    
    #Get the next quantile to remove
    quantile=quantiles[i]
    
    #Create the outfile and set the title for the plot that results from removing the given quantile
    outfile=paste(paste(outdir,paste("corelations_with_3d7_quantile",quantile,sep="_"),sep="/"),"pdf",sep=".")
    
    title=paste("Correlations with 3D7",paste(paste(quantile,"th",sep=""), "Quantile Filtered",sep=" "),sep="\n")
    
    #Get the names of the genes in the given quantile
    genes_in_quantile=names(abs_rank_quantiles)[which(abs_rank_quantiles==quantile)]
    
    #Only retain the data for the genes that are not in the given expression profile
    m_fpkm_data=m_fpkm_data[which(!(rownames(m_fpkm_data) %in% genes_in_quantile)),]
    m_ref_data=m_ref_data[which(!(rownames(m_ref_data) %in% genes_in_quantile)),]

    
    #Get the Correlations of the samples with the reference and put them into a form ggplot can use
    results_cor=cor(x=m_fpkm_data,m_ref_data,method="spearman")
    colnames(results_cor)=seq(1,length(colnames(results_cor)),1)
    df_melted=melt(results_cor)
    
    #Just renaming stuff. This is very specific to the current workflow and should be removed 
    colnames(df_melted)[1]="Strain"
    nf54_string=paste("NF54",timepoint,sep=".")
    pb58_string=paste("PB58",timepoint,sep=".")
    
    #Update how well the 2 samples are correlating with each other
  
    old_ref_cor=cor(m_fpkm_data[,sample2_name],m_fpkm_data[,sample1_name],method="spearman")
    
    #Get the rows corresponding to the timepoints used (the current dataset contains multiple timepoints and
    #only a single timepoint is of interest for any given comparision)
    timepoint_rows=get_rows_coorsponding_to_samples(df_melted,timepoint)
   
    df_timepoint=df_melted[timepoint_rows,]
   
    
    new_first_col=ifelse(df_timepoint[,1]==nf54_string,"K13 WT", "K13 Mutant")
 
    df_timepoint$Strain=new_first_col
    
    
    #Create the correlation plot
    line_plot=ggplot(df_timepoint,aes(x=Var2,y=value,color=Strain))+
      geom_line(size=2)+
      geom_point()+
      ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=26.4,face="bold"))+
      scale_x_continuous(name="3D7 Time Point")+ylab("Spearman Correlation")+
      scale_color_manual(values=c("#ff1aff","#b800e6"))+guides(fill=guide_legend(title="Strain"))
    line_plot=line_plot+theme(axis.text = element_text(size=26),axis.title=element_text(size=26))
    line_plot=line_plot+theme(legend.text=element_text(size=20),legend.key.size=unit(1.5,"cm"))
    pdf(outfile,width=8,height=4.5)
    print(line_plot)
    dev.off()
    
    i=i+1
    
    
  }
  
  
  #Clear the workspace
  rm(list=ls())
  
}

##########Start Analysis

#Infiles

ref_data_infile="Derisi_3d7_ref_data_avg.txt"
sample_fpkm_infile="experimental_data_averages.txt"

#Read in the data
df_ref_data=read.delim(ref_data_infile)
df_fpkm=read.delim(sample_fpkm_infile)

#Format the data for the script
m_ref_data=as.matrix(df_ref_data[,2:length(df_ref_data)])
rownames(m_ref_data)=df_ref_data[,1]
m_fpkm_data=as.matrix(df_fpkm[,2:length(df_fpkm)])
rownames(m_fpkm_data)=df_fpkm[,1]
####Run on real data
 # identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
 #                                     m_ref_data=m_ref_data,
 #                                     quantiles=seq(0,1,0.01),
 #                                     sample1_name="NF54.6h",sample2_name="PB58.6h",
 #                                     minimum_acceptable_correlation=0.8,timepoint="6h",
 #                                     rank_diff_outfile="rank_diff.csv",
 #                                     abs_quantiles_outfile="abs_values_quantiles.csv",
 #                                     rank_diff_hist_outfile="rank_diff_histogram.pdf",
 #                                     outdir="6hr_Correlation_Improvement_Graphs",
 #                                     randomize_data=FALSE,
 #                                     randomization_seed=5712)


####Run on randomized data
identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
                                     m_ref_data=m_ref_data,
                                     quantiles=seq(0,1,0.01),
                                    sample1_name="NF54.6h",sample2_name="PB58.6h",
                                    minimum_acceptable_correlation=0.8,timepoint="6h",
                                    rank_diff_outfile="randomized_rank_diff.csv",
                                    abs_quantiles_outfile="randomized_abs_values_quantiles.csv",
                                    rank_diff_hist_outfile="randomized_rank_diff_histogram.pdf",
                                    outdir="randomized_6hr_Correlation_Improvement_Graphs",
                                    randomize_data=TRUE,
                                    randomization_seed=5712)
