##### description: Gene set expression change trends in M90 and M180-1

#load packages
library(readxl)
library(RColorBrewer)
library(ggplot2)

#Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

#Load count data of 90-day experiment (M90),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)
combat_mat_01 <- read.table("90-day spaceflight_M90/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/RNA-seq/SampleSheet.csv")
combat_mat_01_log <- log2(combat_mat_01+1)
#Load count data of 180-day experiment (M180-1),nobatch_expression_profile.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_180-1, RNA-seq)
combat_mat_02 <- read.table("90-day spaceflight_M180-1/RNA-seq/nobatch_expression_profile.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("90-day spaceflight_M180-1/RNA-seq/SampleSheet.csv")
combat_mat_02_log <- log2(combat_mat_02+1)


pattern_gene_list_01 <- readRDS(file="Differential methylation analysis/result/90-day spaceflight_M90/pattern_gene_list.Rds")
pattern_gene_list_02 <- readRDS(file="Differential methylation analysis/result/180-day spaceflight_M180-1/pattern_gene_list.Rds")


DEG_list_trends_01 <- readRDS("Differential methylation analysis/result//result/DEG_list_trends.Rds")
DEG_list_trends_02 <- readRDS("Differential methylation analysis/result/180-day spaceflight_M180-1/DEG_list_trends.Rds")

library(RColorBrewer)
#term <- "bone resorption"
outdir <- "Multi-group comparisons of geneset expression/result/"
for(term in geneset_all_df_PPI$Description){
   interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
   interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))
   term <- stringr::str_split_i(term, " - ", 1)
   sample_order <- order(SampleSheet_01$Sample_Group,SampleSheet_01$Sample_Well)
   sub_exp_mat <- combat_mat_01_log[interested_geneset,sample_order]
  anno_col <- SampleSheet_01[,c("Sample_Group", "Sample_Well")]
  rownames(anno_col) <-  paste0("90_",SampleSheet_01$Sample_Name)
  exp_trend <- mapply(function(x){names(DEG_list_trends_01[mapply(function(y){x%in%y}, DEG_list_trends_01)])}, interested_geneset)
  exp_trend <- exp_trend[mapply(length,exp_trend)!=0]
  meth_trend <- mapply(function(x){names(pattern_gene_list_01[mapply(function(y){x%in%y}, pattern_gene_list_01)])}, interested_geneset)
  meth_trend <- meth_trend[mapply(length,meth_trend)!=0]
  anno_row <- data.frame(matrix(ncol = 2, nrow = nrow(sub_exp_mat)))
  rownames(anno_row) <- rownames(sub_exp_mat)
  colnames(anno_row) <- c("meth_trend", "exp_trend")
  anno_row$meth_trend[match(names(meth_trend),rownames(anno_row))] <- as.character(meth_trend)
  anno_row$meth_trend[is.na(anno_row$meth_trend)] <- "stable_stable"
  anno_row$exp_trend <- "stable_stable"
  anno_row$exp_trend[match(names(exp_trend),rownames(anno_row))] <- as.character(exp_trend)
  outfile <- paste0("exp_flight01_",term,".pdf")
  my_palette <- c(colorRampPalette(c("#536BE1", "white"))(100),colorRampPalette(c("white", "#FE4C5B"))(100))
  my_color <- brewer.pal(9, "Pastel2")
  gaps_col <- c(3,6)
  num_gene1 <- nrow(sub_exp_mat)
  if(num_gene1<100){
         pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, annotation_colors = list(
           Sample_Group = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Sample_Well = c("H01" = my_color[4], "H02" = my_color[5], "H03" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7])))
  }else{
  pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, width =  1*log(num_gene1), height =  2.5*log2(num_gene1), annotation_colors = list(
           Sample_Group = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Sample_Well = c("H01" = my_color[4], "H02" = my_color[5], "H03" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7])))
    }


  sample_order <- order(SampleSheet_02$Sample_Group,SampleSheet_02$Sample_Well)
  sub_exp_mat <- combat_mat_02_log[interested_geneset,sample_order]
  anno_col <- SampleSheet_02[,c("Sample_Group", "Sample_Well")]
  rownames(anno_col) <- paste0("180_",SampleSheet_02$Sample_Name)
  exp_trend <- mapply(function(x){names(DEG_list_trends_02[mapply(function(y){x%in%y}, DEG_list_trends_02)])}, interested_geneset)
  exp_trend <- exp_trend[mapply(length,exp_trend)!=0]
  meth_trend <- mapply(function(x){names(pattern_gene_list_02[mapply(function(y){x%in%y}, pattern_gene_list_02)])}, interested_geneset)
  meth_trend <- meth_trend[mapply(length,meth_trend)!=0]
  anno_row <- data.frame(matrix(ncol = 2, nrow = nrow(sub_exp_mat)))
  rownames(anno_row) <- rownames(sub_exp_mat)
  colnames(anno_row) <- c("meth_trend", "exp_trend")
  anno_row$meth_trend[match(names(meth_trend),rownames(anno_row))] <- as.character(meth_trend)
  anno_row$meth_trend[is.na(anno_row$meth_trend)] <- "stable_stable"
  anno_row$exp_trend <- "stable_stable"
  anno_row$exp_trend[match(names(exp_trend),rownames(anno_row))] <- as.character(exp_trend)
  outfile <- paste0("exp_flight02_",term,".pdf")
  my_color <- brewer.pal(9, "Pastel2")
  gaps_col <- c(3,6)
  num_gene1 <- nrow(sub_exp_mat)
  if(num_gene1<100){
         pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, annotation_colors = list(
           Sample_Group = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Sample_Well = c("H01" = my_color[4], "H02" = my_color[5], "H03" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7])))
  }else{
  pheatmap::pheatmap(sub_exp_mat, annotation_col = anno_col, annotation_row = anno_row,cluster_cols = F, scale = "row", color = my_palette, filename = paste0(outdir,outfile), gaps_col = gaps_col, width =  1*log(num_gene1), height =  2.5*log2(num_gene1), annotation_colors = list(
           Sample_Group = c("T1" = my_color[1], "T2" = my_color[2], "T3" = my_color[3]),
           Sample_Well = c("H01" = my_color[4], "H02" = my_color[5], "H03" = my_color[6]),
		   meth_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7]),
		   exp_trend = c("down_up" = my_color[1], "stable_stable" = my_color[9], "up_stable" = my_color[3], "down_stable" = my_color[4], "up_down" = my_color[5], "stable_up" = my_color[6], "stable_down" = my_color[2], "up_up" = my_color[8], "down_down" = my_color[7])))
    }

}