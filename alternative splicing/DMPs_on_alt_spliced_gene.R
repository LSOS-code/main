##### description: Statistics on the number of hypermethylated and hypomethylated DMPs on alternatively spliced（AS） genes #####

#Load DMP information
df_DMP_list <- readRDS("Differential methylation analysis/result/DMP_res_M90_M180-1.Rds")
df_DMP_01_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_01_T2_to_T3 <- df_DMP_list[[2]]
df_DMP_01_T1_to_T3 <- df_DMP_list[[3]]
df_DMP_02_T1_to_T2 <- df_DMP_list[[4]]
df_DMP_02_T2_to_T3 <- df_DMP_list[[5]]
df_DMP_02_T1_to_T3 <- df_DMP_list[[6]]

##### 01Start: Obtaining DMPs on alternatively spliced（AS） genes during T1_to_T2 of experiment 01 #####
spling_diff_gene_01_T1_to_T2 <- readxl::read_xlsx(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T2_01_splice_diff.xlsx"))
spling_diff_gene_01_T1_to_T2 <- spling_diff_gene_01_T1_to_T2[,1:16]
splicing_gene_meth_01_T1_to_T2 <- read.table(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T2_01_splice_diff_probes.txt"))
spling_diff_gene_01_T1_to_T2$meth_probes <- splicing_gene_meth_01_T1_to_T2$V3
splicing_gene_meth_list_01_T1_to_T2 <- apply(spling_diff_gene_01_T1_to_T2, 1, function(x){
  gene <- x[[3]]
  probes <- unlist(strsplit(x[[17]],","))[-1]
  probe_anno <- df_DMP_01_T1_to_T2[probes, ]
  probe_anno$gene <- as.character(probe_anno$gene)
  probe_gene <- names(which.max(table(probe_anno$gene)))
  if(sum(probe_anno$gene==gene)>0){
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }else{
    print(c(gene, probe_gene))
    probe_anno$gene[probe_anno$gene==probe_gene] <- gene
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }
  print(table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")]))
  probe_anno$test_id <- x[[1]]
  probe_anno$gene_id <- x[[2]]
  return(probe_anno)
})
lapply(splicing_gene_meth_list_01_T1_to_T2, function(probe_anno){table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05=="up")])})
splicing_gene_meth_df_01_T1_to_T2 <- do.call(rbind, splicing_gene_meth_list_01_T1_to_T2)
#Hypergeometric test for the significance of methylation changes occurring in probes on alternatively spliced genes
N = length(df_DMP_01_T1_to_T2$change_deltaBeta0.05[df_DMP_01_T1_to_T2$change_deltaBeta0.05=="down"])
M = table(df_DMP_01_T1_to_T2$feature[df_DMP_01_T1_to_T2$change_deltaBeta0.05=="down"])["Body"]
n = sum(splicing_gene_meth_df_01_T1_to_T2$change_deltaBeta0.05=="down")
k = sum(splicing_gene_meth_df_01_T1_to_T2$change_deltaBeta0.05=="down" & splicing_gene_meth_df_01_T1_to_T2$feature=="Body")
k/n

q = k
k = n
n = N-M
m = M

p_body_01_T1_to_T2 <- phyper(q-1, m, n, k, lower.tail=F) 
##### 01End: Obtaining DMPs on alternatively spliced（AS） genes during T1_to_T2 of experiment 01 #####


##### 02Start: Obtaining DMPs on alternatively spliced（AS） genes during T2_to_T3 of experiment 01 ######
spling_diff_gene_01_T2_to_T3 <- readxl::read_xlsx(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T2_to_T3_01_splice_diff.xlsx"))
spling_diff_gene_01_T2_to_T3 <- spling_diff_gene_01_T2_to_T3[,1:16]
splicing_gene_meth_01_T2_to_T3 <- read.table(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T2_to_T3_01_splice_diff_probes.txt"))
spling_diff_gene_01_T2_to_T3$meth_probes <- splicing_gene_meth_01_T2_to_T3$V3
splicing_gene_meth_list_01_T2_to_T3 <- apply(spling_diff_gene_01_T2_to_T3, 1, function(x){
  gene <- x[[3]]
  probes <- unlist(strsplit(x[[17]],","))[-1]
  probe_anno <- df_DMP_01_T2_to_T3[probes, ]
  probe_anno$gene <- as.character(probe_anno$gene)
  probe_gene <- names(which.max(table(probe_anno$gene)))
  if(sum(probe_anno$gene==gene)>0){
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }else{
    print(c(gene, probe_gene))
    probe_anno$gene[probe_anno$gene==probe_gene] <- gene
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  } 
  print(table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")]))
  probe_anno$test_id <- x[[1]]
  probe_anno$gene_id <- x[[2]]
  return(probe_anno)
})
lapply(splicing_gene_meth_list_01_T2_to_T3, function(probe_anno){table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")])})

splicing_gene_meth_df_01_T2_to_T3 <- do.call(rbind, splicing_gene_meth_list_01_T2_to_T3)

#Hypergeometric test for the significance of methylation changes occurring in probes on alternatively spliced genes
N = length(df_DMP_01_T2_to_T3$change_deltaBeta0.05[df_DMP_01_T2_to_T3$change_deltaBeta0.05!="stable"])
M = table(df_DMP_01_T2_to_T3$feature[df_DMP_01_T2_to_T3$change_deltaBeta0.05!="stable"])["Body"]
n = sum(splicing_gene_meth_df_01_T2_to_T3$change_deltaBeta0.05!="stable")
k = sum(splicing_gene_meth_df_01_T2_to_T3$change_deltaBeta0.05!="stable" & splicing_gene_meth_df_01_T2_to_T3$feature=="Body")
k/n

q = k
k = n
n = N-M
m = M

p_body_01_T2_to_T3 <- phyper(q-1, m, n, k, lower.tail=F) 
##### 02End: Obtaining DMPs on alternatively spliced（AS） genes for T2_to_T3 of experiment 01 ######


##### 03Start: Obtaining DMPs on alternatively spliced（AS） genes for T1_to_T2 of experiment 02 #######
spling_diff_gene_02_T1_to_T2 <- readxl::read_xlsx(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T2_02_splice_diff.xlsx"))
spling_diff_gene_02_T1_to_T2 <- spling_diff_gene_02_T1_to_T2[,1:16]
splicing_gene_meth_02_T1_to_T2 <- read.table(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T2_02_splice_diff_probes.txt"))
spling_diff_gene_02_T1_to_T2$meth_probes <- splicing_gene_meth_02_T1_to_T2$V3
splicing_gene_meth_list_02_T1_to_T2 <- apply(spling_diff_gene_02_T1_to_T2, 1, function(x){
  gene <- x[[3]]
  probes <- unlist(strsplit(x[[17]],","))[-1]
  probe_anno <- df_DMP_02_T1_to_T2[probes, ]
  probe_anno$gene <- as.character(probe_anno$gene)
  probe_gene <- names(which.max(table(probe_anno$gene)))
  if(sum(probe_anno$gene==gene)>0){
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }else{
    print(c(gene, probe_gene))
    probe_anno$gene[probe_anno$gene==probe_gene] <- gene
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }
  print(table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")]))
  probe_anno$test_id <- x[[1]]
  probe_anno$gene_id <- x[[2]]
  return(probe_anno)
})
lapply(splicing_gene_meth_list_02_T1_to_T2, function(probe_anno){table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")])})
splicing_gene_meth_df_02_T1_to_T2 <- do.call(rbind, splicing_gene_meth_list_02_T1_to_T2)

N = length(df_DMP_02_T1_to_T2$change_deltaBeta0.05[df_DMP_02_T1_to_T2$change_deltaBeta0.05!="stable"])
M = table(df_DMP_02_T1_to_T2$feature[df_DMP_02_T1_to_T2$change_deltaBeta0.05!="stable"])["Body"]
n = sum(splicing_gene_meth_df_02_T1_to_T2$change_deltaBeta0.05!="stable")
k = sum(splicing_gene_meth_df_02_T1_to_T2$change_deltaBeta0.05!="stable" & splicing_gene_meth_df_02_T1_to_T2$feature=="Body")
k/n

q = k
k = n
n = N-M
m = M

p_body_02_T1_to_T2 <- phyper(q-1, m, n, k, lower.tail=F) 
##### 03End: Obtaining DMPs on alternatively spliced（AS） genes for T1_to_T2 of experiment 02 #######


##### 04Start: Obtaining DMPs on alternatively spliced（AS） genes for T2_to_T3 of experiment 02 ########
spling_diff_gene_02_T2_to_T3 <- readxl::read_xlsx(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T2_to_T3_02_splice_diff.xlsx"))
spling_diff_gene_02_T2_to_T3 <- spling_diff_gene_02_T2_to_T3[,1:16]
splicing_gene_meth_02_T2_to_T3 <- read.table(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T2_to_T3_02_splice_diff_probes.txt"))
spling_diff_gene_02_T2_to_T3$meth_probes <- splicing_gene_meth_02_T2_to_T3$V3
splicing_gene_meth_list_02_T2_to_T3 <- apply(spling_diff_gene_02_T2_to_T3, 1, function(x){
  gene <- x[[3]]
  probes <- unlist(strsplit(x[[17]],","))[-1]
  probe_anno <- df_DMP_02_T2_to_T3[probes, ]
  probe_anno$gene <- as.character(probe_anno$gene)
  probe_gene <- names(which.max(table(probe_anno$gene)))
  if(sum(probe_anno$gene==gene)>0){
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }else{
    print(c(gene, probe_gene))
    probe_anno$gene[probe_anno$gene==probe_gene] <- gene
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  } 
  print(table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")]))
  probe_anno$test_id <- x[[1]]
  probe_anno$gene_id <- x[[2]]
  return(probe_anno)
})
lapply(splicing_gene_meth_list_02_T2_to_T3, function(probe_anno){table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")])})

splicing_gene_meth_df_02_T2_to_T3 <- do.call(rbind, splicing_gene_meth_list_02_T2_to_T3)

N = length(df_DMP_02_T2_to_T3$change_deltaBeta0.05[df_DMP_02_T2_to_T3$change_deltaBeta0.05!="stable"])
M = table(df_DMP_02_T2_to_T3$feature[df_DMP_02_T2_to_T3$change_deltaBeta0.05!="stable"])["Body"]
n = sum(splicing_gene_meth_df_02_T2_to_T3$change_deltaBeta0.05!="stable")
k = sum(splicing_gene_meth_df_02_T2_to_T3$change_deltaBeta0.05!="stable" & splicing_gene_meth_df_02_T2_to_T3$feature=="Body")
k/n

q = k
k = n
n = N-M
m = M

p_body_02_T2_to_T3 <- phyper(q-1, m, n, k, lower.tail=F) 
##### 04End: Obtaining DMPs on alternatively spliced（AS） genes for T2_to_T3 of experiment 02 ########


##### 05Start: Obtaining DMPs on alternatively spliced（AS） genes for T1_to_T3 of experiment 02 ########
spling_diff_gene_02_T1_to_T3 <- readxl::read_xlsx(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T3_02_splice_diff.xlsx"))
spling_diff_gene_02_T1_to_T3 <- spling_diff_gene_02_T1_to_T3[,1:16]
splicing_gene_meth_02_T1_to_T3 <- read.table(paste0(spling_diff_gene_list_dir, "alternative splicing/result/T1_to_T3_02_splice_diff_probes.txt"))
spling_diff_gene_02_T1_to_T3$meth_probes <- splicing_gene_meth_02_T1_to_T3$V3
splicing_gene_meth_list_02_T1_to_T3 <- apply(spling_diff_gene_02_T1_to_T3, 1, function(x){
  gene <- x[[3]]
  probes <- unlist(strsplit(x[[17]],","))[-1]
  probe_anno <- df_DMP_02_T1_to_T3[probes, ]
  probe_anno$gene <- as.character(probe_anno$gene)
  probe_gene <- names(which.max(table(probe_anno$gene)))
  if(sum(probe_anno$gene==gene)>0){
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  }else{
    print(c(gene, probe_gene))
    probe_anno$gene[probe_anno$gene==probe_gene] <- gene
    probe_anno <- probe_anno[probe_anno$gene==gene,]
    print(dim(probe_anno))
  } 
  print(table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")]))
  probe_anno$test_id <- x[[1]]
  probe_anno$gene_id <- x[[2]]
  return(probe_anno)
})
lapply(splicing_gene_meth_list_02_T1_to_T3, function(probe_anno){table(probe_anno$feature[which(probe_anno$change_deltaBeta0.05!="stable")])})

splicing_gene_meth_df_02_T1_to_T3 <- do.call(rbind, splicing_gene_meth_list_02_T1_to_T3)

N = length(df_DMP_02_T1_to_T3$change_deltaBeta0.05[df_DMP_02_T1_to_T3$change_deltaBeta0.05!="stable"])
M = table(df_DMP_02_T1_to_T3$feature[df_DMP_02_T1_to_T3$change_deltaBeta0.05!="stable"])["Body"]
n = sum(splicing_gene_meth_df_02_T1_to_T3$change_deltaBeta0.05!="stable")
k = sum(splicing_gene_meth_df_02_T1_to_T3$change_deltaBeta0.05!="stable" & splicing_gene_meth_df_02_T1_to_T3$feature=="Body")
k/n

q = k
k = n
n = N-M
m = M

p_body_02_T1_to_T3 <- phyper(q-1, m, n, k, lower.tail=F) 
##### 05End: Obtaining DMPs on alternatively spliced（AS） genes for T1_to_T3 of experiment 02 ########


#Saving the result tables
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "AS_gene_methProbe_01_T1_to_T2")
addWorksheet(wb, "AS_gene_methProbe_01_T2_to_T3")
addWorksheet(wb, "AS_gene_methProbe_02_T1_to_T2")
addWorksheet(wb, "AS_gene_methProbe_02_T2_to_T3")
addWorksheet(wb, "AS_gene_methProbe_02_T1_to_T3")
writeData(wb, "AS_gene_methProbe_01_T1_to_T2", splicing_gene_meth_df_01_T1_to_T2, rowNames = T)
writeData(wb, "AS_gene_methProbe_01_T2_to_T3", splicing_gene_meth_df_01_T2_to_T3, rowNames = T)

writeData(wb, "AS_gene_methProbe_02_T1_to_T2", splicing_gene_meth_df_02_T1_to_T2, rowNames = T)
writeData(wb, "AS_gene_methProbe_02_T2_to_T3", splicing_gene_meth_df_02_T2_to_T3, rowNames = T)
writeData(wb, "AS_gene_methProbe_02_T1_to_T3", splicing_gene_meth_df_02_T1_to_T3, rowNames = T)
saveWorkbook(wb, "alternative splicing/result/splicing_gene_meth_df.xlsx", overwrite = TRUE)

#Barplot of AS gene number 
splicing_diff_num <- c(length(splicing_gene_meth_list_01_T1_to_T2), length(splicing_gene_meth_list_01_T2_to_T3), 0, length(splicing_gene_meth_list_02_T1_to_T2), length(splicing_gene_meth_list_02_T2_to_T3), length(splicing_gene_meth_list_02_T1_to_T3))
exp_names <- c("90-day flight T1_to_T2", "90-day flight T2_to_T3", "90-day flight 02_T1_to_T3", "180-day flight 02_T1_to_T2", "180-day flight 02_T2_to_T3", "180-day flight 02_T1_to_T3")
splicing_diff_num_df <- data.frame(test_name=factor(exp_names, levels=exp_names), splicing_diff_num=splicing_diff_num, Experiments=c(rep("90-day SZ12", 3), rep("180-day SZ13", 3)))
ggplot(data = splicing_diff_num_df, mapping = aes(x = test_name, y = splicing_diff_num, fill = Experiments)) + 
       geom_bar(stat = 'identity') + 
	   geom_text(aes(label = splicing_diff_num), vjust = -0.5) + 
       scale_fill_manual(values = c("#A71462", "#8695C2"))


#Statistics on the number of sites with methylation changes
t01_T1_to_T2 <- table(splicing_gene_meth_df_01_T1_to_T2$feature[splicing_gene_meth_df_01_T1_to_T2$change_deltaBeta0.05!="stable"])
t01_T2_to_T3 <- table(splicing_gene_meth_df_01_T2_to_T3$feature[splicing_gene_meth_df_01_T2_to_T3$change_deltaBeta0.05!="stable"])

t02_T1_to_T2 <- table(splicing_gene_meth_df_02_T1_to_T2$feature[splicing_gene_meth_df_02_T1_to_T2$change_deltaBeta0.05!="stable"])
t02_T2_to_T3 <- table(splicing_gene_meth_df_02_T2_to_T3$feature[splicing_gene_meth_df_02_T2_to_T3$change_deltaBeta0.05!="stable"])
t02_T1_to_T3 <- table(splicing_gene_meth_df_02_T1_to_T3$feature[splicing_gene_meth_df_02_T1_to_T3$change_deltaBeta0.05!="stable"])

splicing_meth_num_mat <- rbind(t01_T1_to_T2,t01_T2_to_T3,t02_T1_to_T2,t02_T2_to_T3,t02_T1_to_T3)
splicing_meth_num_df <- reshape::melt(splicing_meth_num_mat)
p_values <- round(c(p_body_01_T1_to_T2, p_body_01_T2_to_T3, p_body_02_T1_to_T2,  p_body_02_T2_to_T3, p_body_02_T1_to_T3),3)
p_values[3:4] <- c("8.16e-06", "1.27e-07")
ggplot(data = splicing_meth_num_df, mapping = aes(x = X1, y = value, fill = X2)) + geom_bar(stat = 'identity', position = 'stack') +
  annotate(geom = "text", x = 1, y = 42, label = p_values[1])+
  annotate(geom = "text", x = 2, y = 55, label = p_values[2])+
  annotate(geom = "text", x = 3, y = 90, label = "8.16e-06")+
  annotate(geom = "text", x = 4, y = 90, label = "1.27e-07")+
  annotate(geom = "text", x = 4, y = 90, label = p_values[5])

#Drawing pie charts of the distribution of methylation-altered probes in functional regions of genes
splicing_meth_num_df <- splicing_meth_num_df[order(splicing_meth_num_df$X1, splicing_meth_num_df$X2),]
splicing_meth_num_df$percentage <- unlist(tapply(splicing_meth_num_df$value, factor(splicing_meth_num_df$X1), function(x){x/sum(x)}))*100
splicing_meth_num_df$X1 <- gsub("t01_", "90-days flight ", splicing_meth_num_df$X1)
splicing_meth_num_df$X1 <- gsub("t02_", "180-days flight ", splicing_meth_num_df$X1)
plot <- vector("list", length=length(unique(splicing_meth_num_df$X1)))
names(plot) <- unique(splicing_meth_num_df$X1)
t <- tapply(splicing_meth_num_df$value, factor(splicing_meth_num_df$X1), function(x){sum(x)})
names(p_values) <- unique(splicing_meth_num_df$X1)
colors <- c("#BAE4B3", "#FCDDB8", "#FBF9C1", "#BEDAF6", "#FFB3B3", "#D5EEF1","#FCDAD7", "#C4C6D6")
for(i in unique(splicing_meth_num_df$X1)){
  data <- splicing_meth_num_df[splicing_meth_num_df$X1==i,]
  plot[[i]] <- ggplot(data, aes(x = "", y = percentage, fill = as.factor(X2))) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = colors) +
    coord_polar(theta = "y") +
    labs(fill = "Gene Regions") +
    xlab("percentage") +
    ylab(paste0("N=", t[i], ", P(in Body)=", p_values[i])) +
    ggtitle(i)
}
ggpubr::ggarrange(plot[[1]], plot[[2]], plot[[3]], plot[[4]], nrow = 2, ncol=2, common.legend=T, legend="right")

