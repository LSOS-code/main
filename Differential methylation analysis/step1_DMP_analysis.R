##### description: Differential methylation analysis for 90-day space flight (M90 / 01) and  180-day space flight (M180-1 / 02)
#Due to file size limitations，the complete result tables of differential methylation analysis were stored on https://www.spacelifescience.cn/search/detail?id=24, which was named “Results of differential methylation analysis of DNA methylation probes” 

#load packages
library(ChAMP)
library(ggvenn)
library(ggplot2)
library(openxlsx)


##### 01Start: 90-day differential methylation analysis #####

#Load DNA methylation data of 90-day Experiment (M90),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation)
myCombat_01 <- read.table("90-day spaceflight_M90/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/DNA methylation/SampleSheet.csv")
table(SampleSheet_01$Sample_Group)

#Calculate the significance of differences
DMP_res_01 <- champ.DMP(beta = myCombat_01,
                        pheno = SampleSheet_01$Sample_Group,
                        compare.group = NULL,
                        adjPVal = 1,
                        adjust.method = "BH",
                        arraytype = "EPIC")

df_DMP_01_T1_to_T2 <- DMP_res_01$T1_to_T2
df_DMP_01_T2_to_T3 <- DMP_res_01$T2_to_T3
df_DMP_01_T1_to_T3 <- DMP_res_01$T1_to_T3

#Filtering of probes with a fluctuation range of less than 0.05 
r <- apply(myCombat_01, 1, function(x){max(x)-min(x)})
filter_myCombat_01 <- myCombat_01[r>0.05,]
filter_myCombat_01 <- as.data.frame(t(filter_myCombat_01))
dim(filter_myCombat_01)


#Function of consistency analysis: a consistent direction of change across the three samples
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

#Differential methylation probes(DMPs): probes with a filter fluctuation range greater than 0.05, a p-value for significance of difference less than 0.05, and a consistent direction of change across the three samples
index <- which(SampleSheet_01$Sample_Group%in%c("T1", "T2"))
consist_judge_fun_T1_vs_T2 <- apply(filter_myCombat_01, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_01$Sample_Group[index])})
up_meth_01_T1_to_T2 <- intersect(names(consist_judge_fun_T1_vs_T2)[consist_judge_fun_T1_vs_T2=="Smaller"], 
                                 rownames(df_DMP_01_T1_to_T2)[df_DMP_01_T1_to_T2$P.Value<0.05])
down_meth_01_T1_to_T2 <- intersect(names(consist_judge_fun_T1_vs_T2)[consist_judge_fun_T1_vs_T2=="Greater"], 
                                 rownames(df_DMP_01_T1_to_T2)[df_DMP_01_T1_to_T2$P.Value<0.05])
df_DMP_01_T1_to_T2$change <- "stable"
df_DMP_01_T1_to_T2[up_meth_01_T1_to_T2, "change"] <- "up"
df_DMP_01_T1_to_T2[down_meth_01_T1_to_T2, "change"] <- "down"
table(df_DMP_01_T1_to_T2$change)

index <- which(SampleSheet_01$Sample_Group%in%c("T2", "T3"))
consist_judge_fun_T2_vs_T3 <- apply(filter_myCombat_01, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_01$Sample_Group[index])})
up_meth_01_T2_to_T3 <- intersect(names(consist_judge_fun_T2_vs_T3)[consist_judge_fun_T2_vs_T3=="Smaller"], 
                                 rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$P.Value<0.05])
down_meth_01_T2_to_T3 <- intersect(names(consist_judge_fun_T2_vs_T3)[consist_judge_fun_T2_vs_T3=="Greater"], 
                                   rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$P.Value<0.05])
df_DMP_01_T2_to_T3$change <- "stable"
df_DMP_01_T2_to_T3[up_meth_01_T2_to_T3, "change"] <- "up"
df_DMP_01_T2_to_T3[down_meth_01_T2_to_T3, "change"] <- "down"
table(df_DMP_01_T2_to_T3$change)

index <- which(SampleSheet_01$Sample_Group%in%c("T1", "T3"))
consist_judge_fun_T1_vs_T3 <- apply(filter_myCombat_01, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_01$Sample_Group[index])})
up_meth_01_T1_to_T3 <- intersect(names(consist_judge_fun_T1_vs_T3)[consist_judge_fun_T1_vs_T3=="Smaller"], 
                                 rownames(df_DMP_01_T1_to_T3)[df_DMP_01_T1_to_T3$P.Value<0.05])
down_meth_01_T1_to_T3 <- intersect(names(consist_judge_fun_T1_vs_T3)[consist_judge_fun_T1_vs_T3=="Greater"], 
                                 rownames(df_DMP_01_T1_to_T3)[df_DMP_01_T1_to_T3$P.Value<0.05])
df_DMP_01_T1_to_T3$change <- "stable"
df_DMP_01_T1_to_T3[up_meth_01_T1_to_T3, "change"] <- "up"
df_DMP_01_T1_to_T3[down_meth_01_T1_to_T3, "change"] <- "down"
table(df_DMP_01_T1_to_T3$change)
##### 01End: 90-day differential methylation analysis #####


##### 02Start: M180-1 differential methylation analysis #####

#Load DNA methylation data of 180-day Experiment (M180-1),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, DNA methylation)
myCombat_02 <- read.table("180-day spaceflight_M180-1/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/DNA methylation/SampleSheet.csv")
sum(colnames(myCombat_02)==SampleSheet_02$Sample_Name)

#Calculate the significance of differences
DMP_res_02 <- champ.DMP(beta = myCombat_02,
                     pheno = SampleSheet_02$Sample_Group,
                     compare.group = NULL,
                     adjPVal = 1,
                     adjust.method = "BH",
                     arraytype = "EPIC")

df_DMP_02_T1_to_T2 <- DMP_res_02$T2_to_T1
df_DMP_02_T1_to_T2$deltaBeta <- df_DMP_02_T1_to_T2$deltaBeta*(-1)
DMP_res_02$T1_to_T2 <- df_DMP_02_T1_to_T2
df_DMP_02_T2_to_T3 <- DMP_res_02$T2_to_T3
df_DMP_02_T1_to_T3 <- DMP_res_02$T1_to_T3 

#Filtering of probes with a fluctuation range of less than 0.05
r <- apply(myCombat_02, 1, function(x){max(x)-min(x)})
filter_myCombat_02 <- myCombat_02[r>0.05,]
filter_myCombat_02 <- as.data.frame(t(filter_myCombat_02))
dim(filter_myCombat_02)

#Consistency analysis
#Differential methylation probes(DMPs): probes with a filter fluctuation range greater than 0.05, a p-value for significance of difference less than 0.05, and a consistent direction of change across the three samples
index <- which(SampleSheet_02$Sample_Group%in%c("T1", "T2"))
consist_judge_fun_T1_vs_T2 <- apply(filter_myCombat_02, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_02$Sample_Group[index])})
up_meth_02_T1_to_T2 <- intersect(names(consist_judge_fun_T1_vs_T2)[consist_judge_fun_T1_vs_T2=="Smaller"], 
                                 rownames(df_DMP_02_T1_to_T2)[df_DMP_02_T1_to_T2$P.Value<0.05])
down_meth_02_T1_to_T2 <- intersect(names(consist_judge_fun_T1_vs_T2)[consist_judge_fun_T1_vs_T2=="Greater"], 
                                   rownames(df_DMP_02_T1_to_T2)[df_DMP_02_T1_to_T2$P.Value<0.05])
df_DMP_02_T1_to_T2$change <- "stable"
df_DMP_02_T1_to_T2[up_meth_02_T1_to_T2, "change"] <- "up"
df_DMP_02_T1_to_T2[down_meth_02_T1_to_T2, "change"] <- "down"
table(df_DMP_02_T1_to_T2$change)

index <- which(SampleSheet_02$Sample_Group%in%c("T2", "T3"))
consist_judge_fun_T2_vs_T3 <- apply(filter_myCombat_02, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_02$Sample_Group[index])})
up_meth_02_T2_to_T3 <- intersect(names(consist_judge_fun_T2_vs_T3)[consist_judge_fun_T2_vs_T3=="Smaller"], 
                                 rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$P.Value<0.05])
down_meth_02_T2_to_T3 <- intersect(names(consist_judge_fun_T2_vs_T3)[consist_judge_fun_T2_vs_T3=="Greater"], 
                                   rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$P.Value<0.05])
df_DMP_02_T2_to_T3$change <- "stable"
df_DMP_02_T2_to_T3[up_meth_02_T2_to_T3, "change"] <- "up"
df_DMP_02_T2_to_T3[down_meth_02_T2_to_T3, "change"] <- "down"
table(df_DMP_02_T2_to_T3$change)

index <- which(SampleSheet_02$Sample_Group%in%c("T1", "T3"))
consist_judge_fun_T1_vs_T3 <- apply(filter_myCombat_02, 2, function(x){consist_judge_fun(value=x[index], group=SampleSheet_02$Sample_Group[index])})
up_meth_02_T1_to_T3 <- intersect(names(consist_judge_fun_T1_vs_T3)[consist_judge_fun_T1_vs_T3=="Smaller"], 
                                 rownames(df_DMP_02_T1_to_T3)[df_DMP_02_T1_to_T3$P.Value<0.05])
down_meth_02_T1_to_T3 <- intersect(names(consist_judge_fun_T1_vs_T3)[consist_judge_fun_T1_vs_T3=="Greater"], 
                                 rownames(df_DMP_02_T1_to_T3)[df_DMP_02_T1_to_T3$P.Value<0.05])
df_DMP_02_T1_to_T3$change <- "stable"
df_DMP_02_T1_to_T3[up_meth_02_T1_to_T3, "change"] <- "up"
df_DMP_02_T1_to_T3[down_meth_02_T1_to_T3, "change"] <- "down"
table(df_DMP_02_T1_to_T3$change)
##### 02End: M180-1 differential methylation analysis #####


##### 03End: Statistics of differentially methylated probes (DMPs) at different deltaBeta thresholds #####
df_DMP_list <- c(list(df_DMP_01_T1_to_T2), list(df_DMP_01_T2_to_T3), list(df_DMP_01_T1_to_T3), list(df_DMP_02_T1_to_T2), list(df_DMP_02_T2_to_T3), list(df_DMP_02_T1_to_T3))
names(df_DMP_list) <- c("01_T1_to_T2", "01_T2_to_T3","01_T1_to_T3",  "02_T1_to_T2", "02_T2_to_T3",  "02_T1_to_T3")
deltaBeta_cutoff <- c(0.05, 0.1)
for(i in deltaBeta_cutoff){
   for(j in 1:length(df_DMP_list)){
     DMP_df <- df_DMP_list[[j]]
	 DMP_df$logFC <- DMP_df$deltaBeta
	 new_col <- paste0("change_deltaBeta", i)
	 DMP_df[,new_col] <- DMP_df$change
	 deltaBeta_vec <- DMP_df[DMP_df$change!="stable","deltaBeta"]
	 change_vec <- ifelse(deltaBeta_vec>i, "up", ifelse(deltaBeta_vec<(-i), "down", "stable"))
	 DMP_df[DMP_df$change!="stable",new_col] <- change_vec
	 print(table(DMP_df[,new_col]))
	 df_DMP_list[[j]] <- DMP_df
   }
}


probe_annotation_hg19tohg38 <- read.delim("Differential methylation analysis/genome_conver/probe_annotation_hg19tohg38.txt")
df_DMP_01_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_01_T1_to_T2$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T1_to_T2),"MAPINFO_hg38"]
df_DMP_01_T1_to_T2$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T1_to_T2),"CHR_hg38"]
df_DMP_01_T2_to_T3 <- df_DMP_list[[2]]
df_DMP_01_T2_to_T3$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T2_to_T3),"MAPINFO_hg38"]
df_DMP_01_T2_to_T3$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T2_to_T3),"CHR_hg38"]
df_DMP_01_T1_to_T3 <- df_DMP_list[[3]]
df_DMP_01_T1_to_T3$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T1_to_T3),"MAPINFO_hg38"]
df_DMP_01_T1_to_T3$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_01_T1_to_T3),"CHR_hg38"]

df_DMP_02_T1_to_T2 <- df_DMP_list[[4]]
df_DMP_02_T1_to_T2$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T1_to_T2),"MAPINFO_hg38"]
df_DMP_02_T1_to_T2$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T1_to_T2),"CHR_hg38"]
df_DMP_02_T2_to_T3 <- df_DMP_list[[5]]
df_DMP_02_T2_to_T3$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T2_to_T3),"MAPINFO_hg38"]
df_DMP_02_T2_to_T3$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T2_to_T3),"CHR_hg38"]
df_DMP_02_T1_to_T3 <- df_DMP_list[[6]]
df_DMP_02_T1_to_T3$MAPINFO_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T1_to_T3),"MAPINFO_hg38"]
df_DMP_02_T1_to_T3$CHR_hg38 <- probe_annotation_hg19tohg38[rownames(df_DMP_02_T1_to_T3),"CHR_hg38"]

saveRDS(df_DMP_01_T1_to_T2,file = "Differential methylation analysis/result/df_DMP_01_T1_to_T2.Rds")
saveRDS(df_DMP_01_T2_to_T3,file = "Differential methylation analysis/result/df_DMP_01_T2_to_T3.Rds")
saveRDS(df_DMP_01_T1_to_T3,file = "Differential methylation analysis/result/df_DMP_01_T1_to_T3.Rds")

saveRDS(df_DMP_02_T1_to_T2,file = "Differential methylation analysis/result/df_DMP_02_T1_to_T2.Rds")
saveRDS(df_DMP_02_T2_to_T3,file = "Differential methylation analysis/result/df_DMP_02_T2_to_T3.Rds")
saveRDS(df_DMP_02_T1_to_T3,file = "Differential methylation analysis/result/df_DMP_02_T1_to_T3.Rds")
saveRDS(df_DMP_list, file = "result/DMP_res_M90_M180-1.Rds")

wb <- createWorkbook()
addWorksheet(wb, "DMP_01_T1_to_T2")
addWorksheet(wb, "DMP_01_T2_to_T3")
addWorksheet(wb, "DMP_01_T1_to_T3")
addWorksheet(wb, "DMP_02_T1_to_T2")
addWorksheet(wb, "DMP_02_T2_to_T3")
addWorksheet(wb, "DMP_02_T1_to_T3")
writeData(wb, "DMP_01_T1_to_T2", df_DMP_01_T1_to_T2, rowNames = T)
writeData(wb, "DMP_01_T2_to_T3", df_DMP_01_T2_to_T3, rowNames = T)
writeData(wb, "DMP_01_T1_to_T3", df_DMP_01_T1_to_T3, rowNames = T)
writeData(wb, "DMP_02_T1_to_T2", df_DMP_02_T1_to_T2, rowNames = T)
writeData(wb, "DMP_02_T2_to_T3", df_DMP_02_T2_to_T3, rowNames = T)
writeData(wb, "DMP_02_T1_to_T3", df_DMP_02_T1_to_T3, rowNames = T)
saveWorkbook(wb, "Differential methylation analysis/result/DMP_res_M90_M180-1.xlsx", overwrite = TRUE)


