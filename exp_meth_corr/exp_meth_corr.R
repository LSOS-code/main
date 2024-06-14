##### description: Correlation of LogFCs of gene expression values and promoter methylation levels

#load packages
library(ggplot2)
library(openxlsx)

##### 01Start: Load DMP tables and DEG tables #####
#Load DMP information
df_DMP_list <- readRDS("Differential methylation analysis/result/DMP_res_M90_M180-1.Rds")
df_DMP_01_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_01_T2_to_T3 <- df_DMP_list[[2]]

df_DMP_02_T1_to_T2 <- df_DMP_list[[4]]
df_DMP_02_T2_to_T3 <- df_DMP_list[[5]]

#Load DEG information
DEG_list_trends_01 <- readRDS(file = "DEG_analysis/result/90-day spaceflight_M90/DEG_list_trends.Rds")
DEG_T1_to_T3_01 <- read.table( file = "DEG_analysis/result/90-day spaceflight_M90/T1_to_T3_DEG_table.txt")
DEG_T1_to_T2_01 <- read.table( file = "DEG_analysis/result/90-day spaceflight_M90/T1_to_T2_DEG_table.txt")
DEG_T2_to_T3_01 <- read.table( file = "DEG_analysis/result/90-day spaceflight_M90/T2_to_T3_DEG_table.txt")
DEG_T1_to_T3_01[unlist(DEG_list_trends_01), "change_process"] <-   unlist(lapply(names(DEG_list_trends_01), function(x){rep(x, length(DEG_list_trends_01[[x]]))}))
DEG_T1_to_T3_01$change_process[is.na(DEG_T1_to_T3_01$change_process)] <- "stable_stable"

DEG_list_trends_02 <- readRDS(file = "DEG_analysis/result/180-day spaceflight_M180-1/DEG_list_trends.Rds")
DEG_T1_to_T3_02 <- read.table( file = "DEG_analysis/result/180-day spaceflight_M180-1/T1_to_T3_DEG_table.txt")
DEG_T1_to_T2_02 <- read.table( file = "DEG_analysis/result/180-day spaceflight_M180-1/T1_to_T2_DEG_table.txt")
DEG_T2_to_T3_02 <- read.table( file = "DEG_analysis/result/180-day spaceflight_M180-1/T2_to_T3_DEG_table.txt")
DEG_T1_to_T3_02[unlist(DEG_list_trends_02), "change_process"] <-   unlist(lapply(names(DEG_list_trends_02), function(x){rep(x, length(DEG_list_trends_02[[x]]))}))
DEG_T1_to_T3_02$change_process[is.na(DEG_T1_to_T3_02$change_process)] <- "stable_stable"
##### 01End: Load DMP tables and DEG tables #####


##### 02Start: Classification of DEGs #####
DEG_T1_to_T3_01$baseline <- "Unchanged"
DEG_T1_to_T3_01$baseline[intersect(grep("_down",DEG_T1_to_T3_01$change_process), which(DEG_T1_to_T3_01$logFC<0))] <- "DEG_over"
DEG_T1_to_T3_01$baseline[intersect(grep("_up",DEG_T1_to_T3_01$change_process), which(DEG_T1_to_T3_01$logFC>0))] <- "DEG_over"
DEG_T1_to_T3_01$baseline[DEG_T1_to_T3_01$change_process%in%c("up_stable","down_stable")] <- "DEG_unrecovered"
DEG_T1_to_T3_01$baseline[DEG_T1_to_T3_01$baseline=="Unchanged"&DEG_T1_to_T3_01$change_process!="stable_stable"] <- "Other DEGs"
table(DEG_T1_to_T3_01$baseline)

DEG_T1_to_T3_02$baseline <- "Unchanged"
DEG_T1_to_T3_02$baseline[intersect(grep("_down",DEG_T1_to_T3_02$change_process), which(DEG_T1_to_T3_02$logFC<0))] <- "DEG_over"
DEG_T1_to_T3_02$baseline[intersect(grep("_up",DEG_T1_to_T3_02$change_process), which(DEG_T1_to_T3_02$logFC>0))] <- "DEG_over"
DEG_T1_to_T3_02$baseline[DEG_T1_to_T3_02$change_process%in%c("up_stable","down_stable")] <- "DEG_unrecovered"
DEG_T1_to_T3_02$baseline[DEG_T1_to_T3_02$baseline=="Unchanged"&DEG_T1_to_T3_02$change_process!="stable_stable"] <- "Other DEGs"
table(DEG_T1_to_T3_02$baseline)
##### 02End: Classification of DEGs #####


##### 03Start: Correlation of LogFCs of gene expression values and promoter methylation levels for T2 vs T1 #####
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")
res_df <- c()
#For all genes in M90
all_promoter_probes <- rownames(df_DMP_01_T1_to_T2)[df_DMP_01_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_01))
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "All genes in M90 T2 vs T1")

#For all genes in M180-1
all_promoter_probes <- rownames(df_DMP_02_T1_to_T2)[df_DMP_02_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_02))
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "All genes in M180-1 T2 vs T1")

#For DEGs in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline!="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T1_to_T2)[df_DMP_01_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_01))
genes <- intersect(genes,genes0)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "DEGs in M90 T2 vs T1")

#For DEGs in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline!="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T1_to_T2)[df_DMP_02_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_02))
genes <- intersect(genes,genes0)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "DEGs in M180-1 T2 vs T1")

#For Unchanged genes in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline=="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T1_to_T2)[df_DMP_01_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_01))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Unchanged genes in M90 T2 vs T1")

#For Unchanged genes in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline=="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T1_to_T2)[df_DMP_02_T1_to_T2$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T1_to_T2$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T1_to_T2[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T1_to_T2_02))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T1_to_T2_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Unchanged genes in M180-1 T2 vs T1")
##### 03End: Correlation of LogFCs of gene expression values and promoter methylation levels for T2 vs T1 #####


##### 04Start: Correlation of LogFCs of gene expression values and promoter methylation levels for T3 vs T2 #####
#For DEGs in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline!="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_01))
genes <- intersect(genes,genes0)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "DEGs in M90 T3 vs T2")

#For DEGs in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline!="Unchanged"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_02))
genes <- intersect(genes,genes0)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "DEGs in M180-1 T3 vs T2")

#For over-range rebound DEGs(DEG_over) in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline=="DEG_over"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_01))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Over-range rebound DEGs in M90 T3 vs T2")

#For (DEG_over) in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline=="DEG_over"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_02))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Over-range rebound DEGs in M180-1 T3 vs T2")

#For unrecovered DEGs(DEG_unrecovered) in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline=="DEG_unrecovered"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_01))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Unrecovered DEGs in M90 T3 vs T2")

#For unrecovered DEG(DEG_unrecovered) in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline=="DEG_unrecovered"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_02))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Unrecovered DEGs in M180-1 T3 vs T2")

#For other DEGs in M90
genes0 <- rownames(DEG_T1_to_T3_01)[DEG_T1_to_T3_01$baseline=="Other DEGs"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_01_T2_to_T3)[df_DMP_01_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_01) & df_DMP_01_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_01_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_01))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_01[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Other DEGs in M90 T3 vs T2")

#For other DEGs in M180-1
genes0 <- rownames(DEG_T1_to_T3_02)[DEG_T1_to_T3_02$baseline=="Other DEGs"]
promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")

all_promoter_probes <- rownames(df_DMP_02_T2_to_T3)[df_DMP_02_T2_to_T3$gene%in%rownames(DEG_T1_to_T2_02) & df_DMP_02_T2_to_T3$feature%in%promoter_region]
sub_probe_annotation <- df_DMP_02_T2_to_T3[all_promoter_probes, c("gene", "feature","deltaBeta")]

gene_promoter_deltaBeta <- tapply(sub_probe_annotation$deltaBeta, as.character(sub_probe_annotation$gene), median)
genes <- intersect(names(gene_promoter_deltaBeta), rownames(DEG_T2_to_T3_02))
genes <- intersect(genes,genes0)
length(genes)
gene_promoter_deltaBeta <- gene_promoter_deltaBeta[genes]
gene_expression_logFC <- DEG_T2_to_T3_02[genes,"logFC"]

df <- data.frame(gene_expression_logFC=gene_expression_logFC, gene_promoter_deltaBeta=gene_promoter_deltaBeta)
res <- cor.test(gene_expression_logFC, gene_promoter_deltaBeta)
res_df <- rbind(res_df, c(res$statistic, res$estimate, res$p.value, res$parameter, c(res$conf.int)))
rownames(res_df)[nrow(res_df)] <- c( "Other DEGs in M180-1 T3 vs T2")
##### 04End: Correlation of LogFCs of gene expression values and promoter methylation levels for T3 vs T2 #####


##### 05Start: Display of results #####
#Export correlation analysis results table
colnames(res_df) <- c("statistic","coefficient","P-value","parameter","conf.int_bottom ","conf.int_top")
res_df <- format(res_df, scientific = TRUE, digits = 3)
res_df <- as.data.frame(res_df)
res_df$coefficient <- as.numeric(res_df$coefficient)
res_df$`P-value` <- as.numeric(res_df$`P-value`)
res_df <- as.data.frame(res_df)
res_df$gene_group <- rownames(res_df)
write.csv(res_df, "exp_meth_corr/result/logFC_exp_meth_corr.csv")

#Plotting bar charts to show the results of correlation analyses
p_T2_vs_T1 <- ggplot(res_df[1:6,],aes(x=gsub(" T2 vs T1","",gene_group),y=coefficient,fill=`P-value`))+
    geom_bar(stat = "identity")+geom_text(aes(label = paste0("r = ",coefficient)), vjust = 1, size = 3)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	geom_text(aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3)+
	labs(title="LogFC correlation of gene expression values and promoter methylation levels for T2 vs T1", x="gene group")

p_T3_vs_T2 <- ggplot(res_df[7:14,],aes(x=gsub(" T3 vs T2","",gene_group),y=coefficient,fill=`P-value`))+
    geom_bar(stat = "identity")+geom_text(aes(label = paste0("r = ",coefficient)), vjust = 1, size = 3)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
	geom_text(aes(label = paste0("P = ", `P-value`)), vjust = 2, size = 3)+
	labs(title="LogFC correlation of gene expression values and promoter methylation levels for T3 vs T2", x="gene group")
##### 05End: Display of results #####	
