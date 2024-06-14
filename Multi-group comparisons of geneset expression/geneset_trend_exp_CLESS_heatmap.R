##### description: Gene set expression change trends in 180-day CELSS

#load packages
library(readxl)
library(RColorBrewer)
library(ggplot2)

#Load genesets enriched by PPI nodes 
go_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_GOenrich_res.txt", sep="\t")
kegg_enrich_res <- read.table( file = "WGCNA_PPI_analysis/result/PPI_analysis/function_enrich_res/PPI_node_KEGGenrich_res.txt", sep="\t")
geneset_all_df_PPI <- rbind(go_enrich_res, kegg_enrich_res)

#Load DNA methylation data of 180-days space simulation,nobatch_beta_Mars.txt and SampleSheet.xlsx are available for download at https://www.spacelifescience.cn/search/ (search for 180-days space simulation, DNA methylation)
combat_exp_case <- read.table("180-days space simulation/RNA-seq/nobatch_expression_profile_CELSS.txt", sep="\t", header=T, row.names=1)
sample_info_case <- read.table("180-days space simulation/RNA-seq/samplesheet.txt", sep="\t", header=T)
sample_order <- order(sample_info_case$time)
sorted_time <- c("Pre30", "R2", "R30", "R60", "R75", "R90", "R105", "R120", "R150", "R175", "Post30" )

term = "blood coagulation"
#term = "Osteoclast differentiation"
#term = "NF-kappa B signaling pathway"
#term = "Platelet activation"
interested_geneset <- geneset_all_df_PPI[geneset_all_df_PPI$Description==term, "geneID"]
interested_geneset <- unique(unlist(strsplit(interested_geneset,"/")))
interested_geneset <- interested_geneset[interested_geneset%in%rownames(combat_exp_case)]

PPI_geneset_exp_case <- combat_exp_case[interested_geneset, ]	
annotation_col <- data.frame(time=sample_info_case$time, person=sample_info_case$person)
rownames(annotation_col) <- colnames(PPI_geneset_exp_case)
#pheatmap::pheatmap(PPI_geneset_exp_case, cluster_cols = F, scale = "row", annotation_col = annotation_col, gaps_col = c(1:10*4), cellwidth=10)	
geneset <- rownames(PPI_geneset_exp_case)
times <- unique(sample_info_case$time)
time_num <- length(unique(sample_info_case$time))
annotation <- c()
for (gene in geneset) {
  gene_anno <- c()
  for (j in 1:time_num) {
    test <- wilcox.test(unlist(PPI_geneset_exp_case[gene, 1:4]), unlist(PPI_geneset_exp_case[gene, (j*4-3):(j*4)]))
    P <- round(test$p.value,3)
    gene_anno <- c(gene_anno, P)
  }
  annotation <- rbind(annotation, gene_anno)
}

PPI_geneset_exp_logFC <- apply(PPI_geneset_exp_case, 1, function(x){tapply(x,sample_info_case$time,mean)/mean(x[1:4])})
file <- paste0("Multi-group comparisons of geneset expression/result/", term,"_exp_180_FD_P.pdf")
pheatmap::pheatmap(t(PPI_geneset_exp_logFC), display_numbers = ifelse(annotation<0.05, "*", ""), fontsize=15, cluster_cols = F, height = 6, width = 6, filename= file)	

