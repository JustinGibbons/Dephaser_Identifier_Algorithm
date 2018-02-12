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

v_fpkm_sample_strings=c("NF54.6h","PB58.6h","NF54.24h","PB58.24h","NF54.38h","PB58.38h",
                        "NF54.48h","PB58.48h")



get_columns_coorsponding_to_samples<-function(df_expression,unique_sample_string,ignore.case=TRUE){
  cn=colnames(df_expression)
  return(grep(unique_sample_string,cn,ignore.case = ignore.case))
}

average_columns_by_unique_sample_strings<-function(df,v_sample_strings,ignore.case=TRUE){
  
  averaged_columns=sapply(v_sample_strings,function(g) rowMeans(df[,get_columns_coorsponding_to_samples(df,g,ignore.case=ignore.case)]))
  return(averaged_columns)
}

percentile<-function(input_vector){
  #Returns a vector of the percentile for the values in input_vector
  #There's probably a more efficeint way to do this than a for loop
  percentile_vector=vector(mode="numeric",length=length(input_vector))
  for(i in seq(1:length(input_vector))){
    value=input_vector[[i]]
    number_of_items_less_than_value=sum(I(input_vector<=value))
    percentile_vector[[i]]=number_of_items_less_than_value
    
  }
  percentile_vector=percentile_vector/length(input_vector)
  
  return(percentile_vector)
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
                                     randomize_removal=FALSE,
                                     fold_change_filter=FALSE,
                                     average_expression_filter=FALSE,
                                     raw_fpkm_average_matrix=NULL,
                                     weighted_rank=FALSE,
                                     weighted_rank_constant=2,
                                     randomization_seed=5712
                                     ){
  
#m_fpkm_data is a matrix of the experimental data
#m_ref_data is a matrix of the reference data used to check the life-cycle correlations
#sample1_name and sample2_name are the column names in m_fpkm_data for the samples whose correlations
  #are being computed
#minimum_acceptable_correlation is stopping criteria for the re-phasing algorithm
#quantiles is a vector indicating which quantiles to make
#timepoint is a string that can be used to filter the samples based on the timepoint indicated in the column name
  #and is very specific to this specific dataset
#If fold_change_filter is TRUE you need to supply the FPKM averages for the samples as a matrix
  
#If randomize_data is TRUE will randomize the ranks of the experimental da

  #Create the out directory
  dir.create(outdir)
  
  #Get the genes originally in the data set to figure out which genes were removed
  v_original_genes=rownames(m_fpkm_data)
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
    
  }else if(isTRUE(fold_change_filter)){
    print("Quantiles based on fold change instead of rank difference")
    
    #only include genes used to plot correlations with reference
    raw_fpkm_average_matrix=raw_fpkm_average_matrix[rownames(raw_fpkm_average_matrix) %in% rownames(m_fpkm_data),]
    View(raw_fpkm_average_matrix)
    View(m_fpkm_data)
    
    rank_diff=log2(raw_fpkm_average_matrix[,sample2_name]/raw_fpkm_average_matrix[,sample1_name])
    
    m_fpkm_ranked=rank_matrix_values(m_fpkm_data)
    v_sample1=m_fpkm_ranked[,sample1_name]
    v_sample2=m_fpkm_ranked[,sample2_name]
    
    df_rankdiff=data.frame(gene=rownames(raw_fpkm_average_matrix),rank.diff=rank_diff)
    View(df_rankdiff)
    write.csv(df_rankdiff,rank_diff_outfile)
    abs_rank_diff=abs(rank_diff)
    
    #Plot out the rank differences as a histogram
    pdf(rank_diff_hist_outfile)
    hist(rank_diff,main="Fold Change Distribution", xlab="Fold Change")
    dev.off()
    
  
  }else if(isTRUE(average_expression_filter)){
    #only include genes used to plot correlations with reference
    raw_fpkm_average_matrix=raw_fpkm_average_matrix[rownames(raw_fpkm_average_matrix) %in% rownames(m_fpkm_data),]
    View(raw_fpkm_average_matrix)
    View(m_fpkm_data)
    
    rank_diff=(raw_fpkm_average_matrix[,sample2_name]+raw_fpkm_average_matrix[,sample1_name])/2
    
    m_fpkm_ranked=rank_matrix_values(m_fpkm_data)
    v_sample1=m_fpkm_ranked[,sample1_name]
    v_sample2=m_fpkm_ranked[,sample2_name]
    
    df_rankdiff=data.frame(gene=rownames(raw_fpkm_average_matrix),rank.diff=rank_diff)
    View(df_rankdiff)
    write.csv(df_rankdiff,rank_diff_outfile)
    abs_rank_diff=abs(rank_diff)
    
    #Plot out the rank differences as a histogram
    pdf(rank_diff_hist_outfile)
    hist(rank_diff,main="Fold Change Distribution", xlab="Fold Change")
    dev.off()
    
  }else{
    print("Not Randomizing Data")
    
    m_fpkm_ranked=rank_matrix_values(m_fpkm_data)
    v_sample1=m_fpkm_ranked[,sample1_name]
    v_sample2=m_fpkm_ranked[,sample2_name]
    
    rank_diff=v_sample2-v_sample1
    if(isTRUE(weighted_rank)){
      gene_average=(raw_fpkm_average_matrix[,sample1_name]+raw_fpkm_average_matrix[,sample2_name])/2
      #print(head(gene_average))
      gene_weights=percentile(gene_average)^weighted_rank_constant
      #rank_diff=(rank_diff-min(rank_diff)*(-1))/(max(rank_diff)-min(rank_diff))
      rank_diff=rank_diff/max(abs(rank_diff))
      print(min(rank_diff))
      print(max(rank_diff))
      rank_diff=rank_diff*gene_weights
      
    }
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
  
  #Create a vector to store the changes in correlation
  v_correlations=vector("numeric",length=length(quantiles)+1)
  old_ref_cor=cor(x=v_sample1,y=v_sample2,method="spearman")
  v_correlations[[1]]=old_ref_cor
  
  
  #Loop to perform the filtering
  i=1
  number_of_genes_removed=0
  number_of_genes=nrow(df_fpkm)
  number_of_quantiles=length(quantiles)
  quantile_size=ceiling(number_of_genes/number_of_quantiles)

  #Continue while data correlation is less than minimum correlation and there is still unfiltered data
  while(old_ref_cor<minimum_acceptable_correlation & nrow(m_fpkm_data) >quantile_size){ 
    
    #Get the next quantile to remove
    if(isTRUE(randomize_removal)){
      genes_in_quantile=rownames(m_fpkm_data)[sample(1:length(rownames(m_fpkm_data)),quantile_size)]
      quantile=number_of_quantiles-i+1
      
      
    }else{
      quantile=quantiles[i]
      #Get the names of the genes in the given quantile
      genes_in_quantile=names(abs_rank_quantiles)[which(abs_rank_quantiles==quantile)]
    }
    
    #Create the outfile and set the title for the plot that results from removing the given quantile
    outfile=paste(paste(outdir,paste("corelations_with_3d7_quantile",quantile,sep="_"),sep="/"),"pdf",sep=".")
    
    title=paste("Correlations with 3D7",paste(paste(quantile,"th",sep=""), "Quantile Filtered",sep=" "),sep="\n")
    
    
    
    #Only retain the data for the genes that are not in the given expression profile
    m_fpkm_data=m_fpkm_data[which(!(rownames(m_fpkm_data) %in% genes_in_quantile)),]
    m_ref_data=m_ref_data[which(!(rownames(m_ref_data) %in% genes_in_quantile)),]
    
    number_of_genes_removed=number_of_genes_removed+length(genes_in_quantile)

    
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
   
    v_correlations[[i+1]]=old_ref_cor
    
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
  #Remove the zero correlations from the vector. The zero values are left over from vector initialization.
  v_correlations=v_correlations[v_correlations != 0]
  
  v_genes_remaining=rownames(m_fpkm_data)
  v_genes_removed=setdiff(v_original_genes,v_genes_remaining)
  print(length(v_genes_removed))
  
  return_list=list("Correlations"=v_correlations,
                   "Genes_Removed"=v_genes_removed,
                   "Genes_Not_Removed"=v_genes_remaining)
  
}

rephasing_alg_dot_plot<-function(sampling_correlations_list,outfile_stem){
  #Creates plots to access the performance of the re-phasing algorthm
  
  title="Correlation Improvements\nNonrandom vs Random Data"
  outfile=paste(paste(outfile_stem,"correlation_improvement",sep="_"),"pdf",sep=".")
  
  df_list=vector(mode="list",length=length(sampling_correlations_list))
  
  for(i in seq(1:length(sampling_correlations_list))){
    Correlations=sampling_correlations_list[[i]]
    sample_name=names(sampling_correlations_list)[[i]]
    
    Iterations=seq(1:length(Correlations))
    df=data.frame(Iterations,Correlations)
    df$Data_Type=sample_name
    #print(head(df))
    df_list[[i]]=df
  }
  df_merged=reshape::merge_all(df_list)
  View(df_merged)
    plot=ggplot(df_merged,aes(x=Iterations,y=Correlations,color=Data_Type))
    plot=plot+geom_point(size=2)+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=23,face="bold"))
    plot=plot+scale_x_continuous(name="Algorithm Iteration")+ylab("Spearman Correlation")
    plot=plot+guides(fill=guide_legend(title="Data Type"))
    plot=plot+theme(axis.text = element_text(size=22),axis.title=element_text(size=22))
    plot=plot+theme(legend.text=element_text(size=15),legend.key.size=unit(1.5,"cm"))
    pdf(outfile,width=8,height=4.5)
    print(plot)
    dev.off()
}


# rephasing_alg_dot_plot<-function(non_random_data_correlations, random_data_correlations,outfile_stem){
#   #Creates plots to access the performance of the re-phasing algorthm
# 
#   title="Correlation Improvements\nNonrandom vs Random Data"
# 
#   #Rate of correlation improvement graph
# 
#   outfile=paste(paste(outfile_stem,"correlation_improvement",sep="_"),"pdf",sep=".")
#   #outfile="6hr_correlation_improvement.pdf"
#   # Iterations=seq(1:length(non_random_data_correlations))
#   # Correlations=non_random_data_correlations
#   # df_nonrandom=data.frame(Iterations,Correlations)
#   # df_nonrandom$Data_Type="Nonrandom"
#   #
#   # Iterations=seq(1:length(random_data_correlations))
#   # Correlations=random_data_correlations
#   # df_random=data.frame(Iterations,Correlations)
#   # df_random$Data_Type="Random"
#   #
#   # df_merged=rbind(df_nonrandom,df_random)
#   #
#   # line_plot=ggplot(df_merged,aes(x=Iterations,y=Correlations,color=Data_Type))+
#   #   geom_point(size=2)+ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=26.4,face="bold"))+
#   #   scale_x_continuous(name="Algorithm Iteration")+ylab("Spearman Correlation")+
#   #   scale_color_manual(values=c("#ff1aff","#b800e6"))
#   # line_plot=line_plot+guides(fill=guide_legend(title="Data Type"))
#   # line_plot=line_plot+theme(axis.text = element_text(size=26),axis.title=element_text(size=26))
#   # line_plot=line_plot+theme(legend.text=element_text(size=20),legend.key.size=unit(1.5,"cm"))
#   # pdf(outfile,width=8,height=4.5)
#   # print(line_plot)
#   # dev.off()
# 
# }

rephasing_alg_bar_plot<-function(v_nonrandom_genes_removed,
                                  v_nonrandom_genes_not_removed,
                                  v_random_genes_removed,
                                  v_random_genes_not_removed,outfile_stem){
  #Creates a bar plot of the fraction of total genes removed by the algorithm
  
  #Get the total number of genes in each category
  number_of_nonrandom_genes_removed=length(v_nonrandom_genes_removed)
  total_number_of_nonrandom_genes=length(union(v_nonrandom_genes_removed,v_nonrandom_genes_not_removed))
  percent_nonrandom_removed=(number_of_nonrandom_genes_removed/total_number_of_nonrandom_genes)*100
  print(percent_nonrandom_removed)
  
  number_of_random_genes_removed=length(v_random_genes_removed)
  total_number_of_random_genes=length(union(v_random_genes_removed,v_random_genes_not_removed))
  percent_random_removed=(number_of_random_genes_removed/total_number_of_random_genes)*100
  print(percent_random_removed)
  
  if(total_number_of_nonrandom_genes != total_number_of_random_genes){
    print("Warning: Total number of genes was not equal on the nonrandom and random groups")
    print("Your graph will not be made!")
  }else{
    outfile=paste(paste(outfile_stem,"barplot",sep="_"),"pdf",sep=".")
    df=data.frame(Data_Type=c("Nonrandom","Random"),
                  Percent_Removed=c(percent_nonrandom_removed,percent_random_removed))

    plot=ggplot(data=df, aes(x=Data_Type,y=Percent_Removed))
    plot=plot+geom_bar(stat="identity")
    
    pdf(outfile)
    print(plot)
    dev.off()
    
  }
  
  
}

volcano_plot<-function(m_average_expression, group1_identifier, outfile_stem=NULL,group2_identifier,flagged_genes,use_rankdiff=FALSE){
  #Creates a volcano plot using the expression data in the supplied matrix
  #Input is a matrix were the rows are genes and the columns are samples. Assumes 1 has already been added
  #to each sample and replicates are averaged.
  #The different sample types which are identied by the variables group1_identifier and group2_identifier.
  #flagged_genes is a vector of the genes that are too be flagged with a different color
  #The fold change is than calculated as log2(group1/group2)
  #xaxis is the average expression in the 2 groups
  
  
  
  if(isTRUE(use_rankdiff)){
    outfile=paste(paste(outfile_stem,group1_identifier,group2_identifier,"rank_diff_volano_plot",sep="_"),"pdf",sep=".")
    ranked_matrix=rank_matrix_values(m_average_expression)
    rank_diff=ranked_matrix[,group1_identifier]-ranked_matrix[,group2_identifier]
    gene_flagged=rownames(m_average_expression) %in% flagged_genes
    average_expression=(m_average_expression[,group1_identifier]+m_average_expression[,group2_identifier])/2
    
    df=data.frame(GeneID=rownames(m_average_expression),Rank_Diff=rank_diff,
                  Average_Expression=log10(average_expression),
                  Gene_Flagged=gene_flagged)
    plot=ggplot(df,aes(x=Average_Expression,y=Rank_Diff,color=Gene_Flagged))
    
    
  }else{
  
  outfile=paste(paste(outfile_stem,group1_identifier,group2_identifier,"fold_change_volano_plot",sep="_"),"pdf",sep=".")
  fold_change=log2(m_average_expression[,group1_identifier]/m_average_expression[,group2_identifier])
  gene_flagged=rownames(m_average_expression) %in% flagged_genes
  average_expression=(m_average_expression[,group1_identifier]+m_average_expression[,group2_identifier])/2
  df=data.frame(GeneID=rownames(m_average_expression),Fold_Change=fold_change,
                Average_Expression=log10(average_expression),
                Gene_Flagged=gene_flagged)
  plot=ggplot(df,aes(x=Average_Expression,y=Fold_Change,color=Gene_Flagged))
  #plot=plot+scale_y_continuous(limits=c(-5,5))
  }
  
  
  
  plot=plot+geom_point(size=2)
  #ggtitle(title)+theme(plot.title = element_text(hjust = 0.5,size=26.4,face="bold"))+
  #   scale_x_continuous(name="Algorithm Iteration")+ylab("Spearman Correlation")+
  #   scale_color_manual(values=c("#ff1aff","#b800e6"))
  # line_plot=line_plot+guides(fill=guide_legend(title="Data Type"))
  # line_plot=line_plot+theme(axis.text = element_text(size=26),axis.title=element_text(size=26))
  # line_plot=line_plot+theme(legend.text=element_text(size=20),legend.key.size=unit(1.5,"cm"))
  pdf(outfile)
  print(plot)
  dev.off()

  
}


##########Start Analysis

#Infiles

ref_data_infile="Derisi_3d7_ref_data_avg.txt"
sample_fpkm_infile="experimental_data_averages.txt"
sample_raw_fpkm_infile="nf54_pb58_gene_fpkm.tab"
#Read in the data
df_ref_data=read.delim(ref_data_infile)
df_fpkm=read.delim(sample_fpkm_infile)
df_raw_fpkm=read.delim(sample_raw_fpkm_infile)

#Format the data for the script
m_ref_data=as.matrix(df_ref_data[,2:length(df_ref_data)])
rownames(m_ref_data)=df_ref_data[,1]
m_fpkm_data=as.matrix(df_fpkm[,2:length(df_fpkm)])
rownames(m_fpkm_data)=df_fpkm[,1]
m_raw_fpkm=as.matrix(df_raw_fpkm[,2:length(df_raw_fpkm)])+1 #Adding 1 to everything to pervent zero division
rownames(m_raw_fpkm)=df_raw_fpkm[,1]
m_raw_filtered_fpkm=m_raw_fpkm[rownames(m_raw_fpkm) %in% rownames(m_fpkm_data),] #Only keep genes that are in m_fpkm_data

#Average the experimental data
v_fpkm_sample_strings=c("NF54.6h","PB58.6h","NF54.24h","PB58.24h","NF54.38h","PB58.38h",
                        "NF54.48h","PB58.48h")
m_raw_filtered_fpkm_avgs=average_columns_by_unique_sample_strings(m_raw_filtered_fpkm,v_fpkm_sample_strings)
m_raw_fpkm_avgs=average_columns_by_unique_sample_strings(m_raw_fpkm,v_fpkm_sample_strings)




####Run on real data
# l_nonrandom_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                      m_ref_data=m_ref_data,
#                                      quantiles=seq(0,1,0.01),
#                                      sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                      minimum_acceptable_correlation=0.8,timepoint="6h",
#                                      rank_diff_outfile="rank_diff.csv",
#                                      abs_quantiles_outfile="abs_values_quantiles.csv",
#                                      rank_diff_hist_outfile="rank_diff_histogram.pdf",
#                                      outdir="6hr_Correlation_Improvement_Graphs",
#                                      randomize_data=FALSE,
#                                      randomization_seed=5712)
nonrandom_correlations=l_nonrandom_results$Correlations
nonrandom_genes_removed=l_nonrandom_results$Genes_Removed
nonrandom_genes_not_removed=l_nonrandom_results$Genes_Not_Removed

# print(head(nonrandom_correlations))
# print(head(nonrandom_genes_removed))
# print(head(nonrandom_genes_not_removed))




####Run on randomized data
# l_randomized_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                      m_ref_data=m_ref_data,
#                                      quantiles=seq(0,1,0.01),
#                                     sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                     minimum_acceptable_correlation=0.8,timepoint="6h",
#                                     rank_diff_outfile="randomized_rank_diff.csv",
#                                     abs_quantiles_outfile="randomized_abs_values_quantiles.csv",
#                                     rank_diff_hist_outfile="randomized_rank_diff_histogram.pdf",
#                                     outdir="randomized_6hr_Correlation_Improvement_Graphs",
#                                     randomize_data=TRUE,
#                                     randomization_seed=5712)

randomized_correlations=l_randomized_results$Correlations
randomized_genes_removed=l_randomized_results$Genes_Removed
randomized_genes_not_removed=l_randomized_results$Genes_Not_Removed

# l_randomized_removal_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                      m_ref_data=m_ref_data,
#                                      quantiles=seq(0,1,0.01),
#                                     sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                     minimum_acceptable_correlation=0.8,timepoint="6h",
#                                     rank_diff_outfile="randomized_rank_diff.csv",
#                                     abs_quantiles_outfile="randomized_removal_abs_values_quantiles.csv",
#                                     rank_diff_hist_outfile="randomized_removal_rank_diff_histogram.pdf",
#                                     outdir="randomized_removal_6hr_Correlation_Improvement_Graphs",
#                                     randomize_data=FALSE,
#                                     randomize_removal=TRUE,
#                                     randomization_seed=5712)

randomized_removal_correlations=l_randomized_removal_results$Correlations
randomized_removal_genes_removed=l_randomized_removal_results$Genes_Removed
randomized_removal_genes_not_removed=l_randomized_removal_results$Genes_Not_Removed

# l_fold_change_removal_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                                                   m_ref_data=m_ref_data,
#                                                                   quantiles=seq(0,1,0.01),
#                                                                   sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                                                   minimum_acceptable_correlation=0.8,timepoint="6h",
#                                                                   rank_diff_outfile="fold_change_rank_diff.csv",
#                                                                   abs_quantiles_outfile="fold_change_removal_abs_values_quantiles.csv",
#                                                                   rank_diff_hist_outfile="fold_change_removal_rank_diff_histogram.pdf",
#                                                                   outdir="fold_change_removal_6hr_Correlation_Improvement_Graphs",
#                                                                   randomize_data=FALSE,
#                                                                   randomize_removal=FALSE,
#                                                                   fold_change_filter = TRUE,
#                                                                   raw_fpkm_average_matrix=m_raw_filtered_fpkm_avgs,
#                                                                   randomization_seed=5712)

fold_change_removal_correlations=l_fold_change_removal_results$Correlations
fold_change_removal_genes_removed=l_fold_change_removal_results$Genes_Removed
fold_change_removal_genes_not_removed=l_fold_change_removal_results$Genes_Not_Removed

# l_average_expression_removal_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                                                    m_ref_data=m_ref_data,
#                                                                    quantiles=seq(0,1,0.01),
#                                                                    sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                                                    minimum_acceptable_correlation=0.8,timepoint="6h",
#                                                                    rank_diff_outfile="average_expression_rank_diff.csv",
#                                                                    abs_quantiles_outfile="average_expression_removal_abs_values_quantiles.csv",
#                                                                    rank_diff_hist_outfile="average_expression_removal_rank_diff_histogram.pdf",
#                                                                    outdir="average_expression_removal_6hr_Correlation_Improvement_Graphs",
#                                                                    randomize_data=FALSE,
#                                                                    randomize_removal=FALSE,
#                                                                    fold_change_filter = FALSE,
#                                                                    average_expression_filter = TRUE,
#                                                                    raw_fpkm_average_matrix=m_raw_filtered_fpkm_avgs,
#                                                                    randomization_seed=5712)

average_expression_removal_correlations=l_average_expression_removal_results$Correlations
average_expression_removal_genes_removed=l_average_expression_removal_results$Genes_Removed
average_expression_removal_genes_not_removed=l_average_expression_removal_results$Genes_Not_Removed

# print(head(randomized_correlations))
# print(head(randomized_genes_removed))
# print(head(randomized_genes_not_removed))

# l_weighted_rank2_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                      m_ref_data=m_ref_data,
#                                      quantiles=seq(0,1,0.01),
#                                      sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                      minimum_acceptable_correlation=0.8,timepoint="6h",
#                                      rank_diff_outfile="weighted_rank2_rank_diff.csv",
#                                      abs_quantiles_outfile="weighted_rank2_abs_values_quantiles.csv",
#                                      rank_diff_hist_outfile="weighted_rank2_rank_diff_histogram.pdf",
#                                      outdir="weighted_rank2_6hr_Correlation_Improvement_Graphs",
#                                      randomize_data=FALSE,
#                                      weighted_rank = TRUE,
#                                      raw_fpkm_average_matrix=m_raw_filtered_fpkm_avgs,
#                                      weighted_rank_constant = 2,
#                                      randomization_seed=5712)

weighted_rank2_correlations=l_weighted_rank2_results$Correlations
weighted_rank2_genes_removed=l_weighted_rank2_results$Genes_Removed
weighted_rank2_genes_not_removed=l_weighted_rank2_results$Genes_Not_Removed

# l_weighted_rank1.5_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                                               m_ref_data=m_ref_data,
#                                                               quantiles=seq(0,1,0.01),
#                                                               sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                                               minimum_acceptable_correlation=0.8,timepoint="6h",
#                                                               rank_diff_outfile="weighted_rank1.5_rank_diff.csv",
#                                                               abs_quantiles_outfile="weighted_rank_abs_values_quantiles.csv",
#                                                               rank_diff_hist_outfile="weighted_rank1.5_rank_diff_histogram.pdf",
#                                                               outdir="weighted_rank1.5_6hr_Correlation_Improvement_Graphs",
#                                                               randomize_data=FALSE,
#                                                               weighted_rank = TRUE,
#                                                               raw_fpkm_average_matrix=m_raw_filtered_fpkm_avgs,
#                                                               weighted_rank_constant = 1.5,
#                                                               randomization_seed=5712)

weighted_rank1.5_correlations=l_weighted_rank1.5_results$Correlations
weighted_rank1.5_genes_removed=l_weighted_rank1.5_results$Genes_Removed
weighted_rank1.5_genes_not_removed=l_weighted_rank1.5_results$Genes_Not_Removed

# l_weighted_rank1_results=identify_dephasing_drivers_nf54_pb58(m_fpkm_data=m_fpkm_data,
#                                                               m_ref_data=m_ref_data,
#                                                               quantiles=seq(0,1,0.01),
#                                                               sample1_name="NF54.6h",sample2_name="PB58.6h",
#                                                               minimum_acceptable_correlation=0.8,timepoint="6h",
#                                                               rank_diff_outfile="weighted_rank1_rank_diff.csv",
#                                                               abs_quantiles_outfile="weighted_rank_abs_values_quantiles.csv",
#                                                               rank_diff_hist_outfile="weighted_rank1_rank_diff_histogram.pdf",
#                                                               outdir="weighted_rank1_6hr_Correlation_Improvement_Graphs",
#                                                               randomize_data=FALSE,
#                                                               weighted_rank = TRUE,
#                                                               raw_fpkm_average_matrix=m_raw_filtered_fpkm_avgs,
#                                                               weighted_rank_constant = 1,
#                                                               randomization_seed=5712)

weighted_rank1_correlations=l_weighted_rank1_results$Correlations
weighted_rank1_genes_removed=l_weighted_rank1_results$Genes_Removed
weighted_rank1_genes_not_removed=l_weighted_rank1_results$Genes_Not_Removed


# l_removing_type=list("Normalized Expression\nRank Difference"=nonrandom_correlations,
#                      "Random Removal"=randomized_removal_correlations,
#                      "Fold Change Removal"=fold_change_removal_correlations,
#                      "Average Expression Removal"=average_expression_removal_correlations)

l_removing_type=list("Normalized Expression\nRank Difference"=nonrandom_correlations,
                     "Random Removal"=randomized_removal_correlations,
                     "Average Expression Removal"=average_expression_removal_correlations)


rephasing_alg_dot_plot(l_removing_type,"6hr_average_expression_removal")

#rephasing_alg_bar_plot(nonrandom_genes_removed,nonrandom_genes_not_removed,
#                       randomized_genes_removed, randomized_genes_not_removed,outfile_stem = "6hr")


#volcano_plot(m_raw_fpkm_avgs,group1_identifier ="PB58.6h" ,group2_identifier = "NF54.6h",flagged_genes=nonrandom_genes_removed) #Want the reference to be
                                                                                          #the denominator

# volcano_plot(m_raw_fpkm_avgs,group1_identifier ="PB58.6h" ,group2_identifier = "NF54.6h",nonrandom_genes_removed,use_rankdiff = TRUE) #Want the reference to be
#the denominator

#volcano_plot(m_raw_fpkm_avgs,group1_identifier ="PB58.6h" ,group2_identifier = "NF54.6h",
#             outfile_stem="percentile_weighted2",weighted_rank2_genes_removed) #Want the reference to be
                                                                                           #the denominator

# volcano_plot(m_raw_fpkm_avgs,group1_identifier ="PB58.6h" ,group2_identifier = "NF54.6h",
#             outfile_stem="percentile_weighted1",weighted_rank1_genes_removed) #Want the reference to be
#                                                                                 #the denominator
# volcano_plot(m_raw_fpkm_avgs,group1_identifier ="PB58.6h" ,group2_identifier = "NF54.6h",
#             outfile_stem="percentile_weighted1.5",weighted_rank1.5_genes_removed) #Want the reference to be
#                                                                                 #the denominator