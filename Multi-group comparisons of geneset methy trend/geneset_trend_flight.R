##### description: Gene set methylation change trends in M90 and M180-1

#load packages
library(readxl)
library(RColorBrewer)
library(ggplot2)


#Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

##############################################################
#Gene set methylation change trends in M90(01)
df_DMP_01_T1_to_T2 <- readRDS("Differential methylation analysis/result/90-day spaceflight_M90/df_DMP_01_T1_to_T2.Rds")
#Load DNA methylation data of 90-day Experiment (M90),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation)
myCombat01 <- read.table("90-day spaceflight_M90/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/DNA methylation/SampleSheet.csv")
plot_list_flight01 <- c()
for(term in geneset_all_df_PPI$Description){
  interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
  interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))
  
  promoter_region <- c("TSS1500","TSS200","5'UTR","1stExon")
  interested_geneset_annotation_01 <- df_DMP_01_T1_to_T2[df_DMP_01_T1_to_T2$gene%in%interested_geneset & df_DMP_01_T1_to_T2$feature%in%promoter_region, ]
  interested_geneset_meth_probe_01 <-  as.data.frame(myCombat01[rownames(interested_geneset_annotation_01), ])
  interested_geneset_meth_probe_01$gene <- interested_geneset_annotation_01$gene

  interested_geneset_meth_mat_01 <- aggregate(interested_geneset_meth_probe_01[,-ncol(interested_geneset_meth_probe_01)], by = list(group = interested_geneset_annotation_01$gene), median)
  rownames(interested_geneset_meth_mat_01) <- interested_geneset_meth_mat_01$group
  interested_geneset_meth_mat_01 <- interested_geneset_meth_mat_01[,-1]

  consist_res_T1T2 <- apply(interested_geneset_meth_mat_01[,SampleSheet_01$Sample_Group%in%c("T1","T2")], 1, consist_judge_fun, group=SampleSheet_01$Sample_Group[SampleSheet_01$Sample_Group%in%c("T1","T2")])
  consist_res_T2T3 <- apply(interested_geneset_meth_mat_01[,SampleSheet_01$Sample_Group%in%c("T2","T3")], 1, consist_judge_fun, group=SampleSheet_01$Sample_Group[SampleSheet_01$Sample_Group%in%c("T2","T3")])

  consist_index <- which(consist_res_T1T2!="Inconsistent"&consist_res_T2T3!="Inconsistent")
  filter_interested_geneset_meth_mat_01 <- interested_geneset_meth_mat_01[consist_index, ]

  if(nrow(filter_interested_geneset_meth_mat_01)>1){
     r <- apply(filter_interested_geneset_meth_mat_01, 1, function(x){max(x)-min(x)})
     filter_interested_geneset_meth_mat_01 <- filter_interested_geneset_meth_mat_01[r>0.05,]
     #print(dim(filter_interested_geneset_meth_mat_01))

     if(nrow(filter_interested_geneset_meth_mat_01)>1){
        mean_meth_mat_01 <- apply(filter_interested_geneset_meth_mat_01, 1, function(x) tapply(x, SampleSheet_01$Sample_Group, mean))
        mean_meth_mat_01 <- t(mean_meth_mat_01)
		print(dim(mean_meth_mat_01))
        trend_p <- trend_plot(mean_meth_mat_01)
        plot_list_flight01 <- c(plot_list_flight01, list(trend_p))
		term <- stringr::str_split_i(term, " - ", 1)
        names(plot_list_flight01)[length(plot_list_flight01)] <- term
     }

    }
}


saveRDS(plot_list_flight01, "Multi-group comparisons of geneset methy trend/result/plot_list_flight01.Rds")
for(i in 1:length(plot_list_flight01)){
    filename<-paste("Multi-group comparisons of geneset methy trend/result/flight01_",names(plot_list_flight01)[i],".pdf",sep = "")
    pdf(filename,width=5,height=3)
	trend_p <- plot_list_flight01[i]
    print(trend_p)
    dev.off()
}


###################################################################################
#Gene set methylation change trends in M180-1(02)
setwd(mydir)
df_DMP_02_T1_to_T2 <- readRDS("Differential methylation analysis/result/90-day spaceflight_M180-1/df_DMP_02_T1_to_T2.Rds")
#Load DNA methylation data of 180-day Experiment (M180-1),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation)
myCombat02 <- read.table("180-day spaceflight_M180-1/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/DNA methylation/SampleSheet.csv")
plot_list_flight02 <- c()
for(term in geneset_all_df_PPI$Description){
  interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
  interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))

  promoter_region <- c("TSS1500","TSS200","5'UTR","1stExon")
  interested_geneset_annotation_02 <- df_DMP_02_T1_to_T2[df_DMP_02_T1_to_T2$gene%in%interested_geneset & df_DMP_02_T1_to_T2$feature%in%promoter_region, ]
  interested_geneset_meth_probe_02 <-  as.data.frame(myCombat02[rownames(interested_geneset_annotation_02), ])
  interested_geneset_meth_probe_02$gene <- interested_geneset_annotation_02$gene

  interested_geneset_meth_mat_02 <- aggregate(interested_geneset_meth_probe_02[,-ncol(interested_geneset_meth_probe_02)], by = list(group = interested_geneset_annotation_02$gene), median)
  rownames(interested_geneset_meth_mat_02) <- interested_geneset_meth_mat_02$group
  interested_geneset_meth_mat_02 <- interested_geneset_meth_mat_02[,-1]

  consist_res_T1T2 <- apply(interested_geneset_meth_mat_02[,SampleSheet_02$Sample_Group%in%c("T1","T2")], 1, consist_judge_fun, group=SampleSheet_02$Sample_Group[SampleSheet_02$Sample_Group%in%c("T1","T2")])
  consist_res_T2T3 <- apply(interested_geneset_meth_mat_02[,SampleSheet_02$Sample_Group%in%c("T2","T3")], 1, consist_judge_fun, group=SampleSheet_02$Sample_Group[SampleSheet_02$Sample_Group%in%c("T2","T3")])

  consist_index <- which(consist_res_T1T2!="Inconsistent"&consist_res_T2T3!="Inconsistent")
  filter_interested_geneset_meth_mat_02 <- interested_geneset_meth_mat_02[consist_index, ]

  if(nrow(filter_interested_geneset_meth_mat_02)>1){
     r <- apply(filter_interested_geneset_meth_mat_02, 1, function(x){max(x)-min(x)})
     filter_interested_geneset_meth_mat_02 <- filter_interested_geneset_meth_mat_02[r>0.05,]
     dim(filter_interested_geneset_meth_mat_02)

     if(nrow(filter_interested_geneset_meth_mat_02)>1){
        mean_meth_mat_02 <- apply(filter_interested_geneset_meth_mat_02, 1, function(x) tapply(x, SampleSheet_02$Sample_Group, mean))
        mean_meth_mat_02 <- t(mean_meth_mat_02)
        print(dim(mean_meth_mat_02))
        trend_p <- trend_plot(mean_meth_mat_02)

        plot_list_flight02 <- c(plot_list_flight02, list(trend_p))
		term <- stringr::str_split_i(term, " - ", 1)
        names(plot_list_flight02)[length(plot_list_flight02)] <- term
     }

    }
}

saveRDS(plot_list_flight02, "Multi-group comparisons of geneset methy trend/result/plot_list_flight02.Rds")

for(i in 1:length(plot_list_flight02)){
    filename<-paste("Multi-group comparisons of geneset methy trend/result/flight02_",names(plot_list_flight02)[i],".pdf",sep = "")
    pdf(filename,width=5,height=3)
	trend_p <- plot_list_flight02[i]
    print(trend_p)
    dev.off()
}
################################

####################################


consist_judge_fun <- function(value, group) {
  value_by_group <- split(value, group)
  #print(paste0("compare between ", paste0(names(value_by_group), collapse = ", ")))
  greater_judge <- value_by_group[[1]] > value_by_group[[2]]
  smaller_judge <- value_by_group[[1]] < value_by_group[[2]]
  if(sum(greater_judge)==3){
    return("Greater")
  }else if(sum(smaller_judge)==3){
    return("Smaller")
  }else{
      return("Inconsistent")
  }
}



#trend_matrix: col:gene; rownames:time by order; value:mean of methlation beta value of the genes.
trend_plot <- function(trend_matrix){
  time<-colnames(trend_matrix)
  trend_matrix<-t(apply(trend_matrix, 1, function(x)scale(x, center = F)))

  t<-1:length(time)
  data<-as.data.frame(t(trend_matrix))
  data$t<-t
  data_long<-reshape::melt(data,id.vars = "t")
  x<-1:length(time)
  quantile_value<-apply(trend_matrix,2,function(x){
    temp<-round(quantile(x,seq(0,1,0.01)),3)
    return(temp)
  })
  pdata.list<-list()
  for(k in 1:50){
    pdata<-data.frame(x,lower = quantile_value[k,],upper = quantile_value[102-k,])
    pdata.list[[k]]<-pdata
  }
  my_col <- c("#C6E4F5", "#2BAAF3")
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend<-ggplot()+
    geom_rect(aes(xmin = -Inf,xmax = 1,ymin = -Inf,ymax = Inf),fill = "#FFFFCC",alpha = 0.4)+
    geom_rect(aes(xmin = 1,xmax = length(time),ymin = -Inf,ymax = Inf),fill = "#DEF6F3",alpha = 0.4)+
    geom_rect(aes(xmin = length(time),xmax = Inf,ymin = -Inf,ymax = Inf),fill = "#FFFFCC",alpha = 0.4)+
    theme(panel.grid=element_blank(),panel.border=element_rect(fill=NA,color="black", size=1))

  for(k in 6:50){
    pdata<-pdata.list[[k]]
    plot.trend<-plot.trend+geom_ribbon(data = pdata,aes(ymin=lower, ymax=upper, x=x), fill = myPalette[k-5], alpha = 1)
  }
  plot.trend<-plot.trend+
    theme_bw()+
    xlab("")+
    theme(panel.grid =element_blank()) + 
    scale_x_continuous(breaks=seq(1, length(time), 1), labels = time)
  return(plot.trend)
}


