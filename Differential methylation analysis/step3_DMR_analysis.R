setwd("F:/database")
myCombat_01 <- read.table("90-day spaceflight_M90/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_01 <- read.csv("90-day spaceflight_M90/DNA methylation/SampleSheet.csv")
table(SampleSheet_01$Sample_Group)

myCombat_02 <- read.table("180-day spaceflight_M180-1/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_02 <- read.csv("180-day spaceflight_M180-1/DNA methylation/SampleSheet.csv")

myCombat_03 <- read.table("180-day spaceflight_M180-2/DNA methylation/nobatch_beta.txt", header=T, row.names=1, as.is=T, sep="\t")
SampleSheet_03 <- read.csv("180-day spaceflight_M180-2/DNA methylation/SampleSheet.csv")


df_DMP_list <- readRDS("Differential methylation analysis/result/DMP_res_M90_M180-1.Rds")
df_DMP_01_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_01_T2_to_T3 <- df_DMP_list[[2]]
df_DMP_01_T1_to_T3 <- df_DMP_list[[3]]

df_DMP_02_T1_to_T2 <- df_DMP_list[[4]]
df_DMP_02_T2_to_T3 <- df_DMP_list[[5]]
df_DMP_02_T1_to_T3 <- df_DMP_list[[6]]
rm(df_DMP_list)

df_DMP_list <- readRDS("Differential methylation analysis/result/DMP_res_M180-2.Rds")
df_DMP_03_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_03_T2_to_T3 <- df_DMP_list[[3]]
df_DMP_03_T1_to_T3 <- df_DMP_list[[2]]
rm(df_DMP_list)


mydir <- "D:/PhD/major_program/analysis/meth/stp2_DMP/DMR/data"
setwd(mydir)
   
pheno_01_T1_to_T2 <- factor(SampleSheet_01$Sample_Group[SampleSheet_01$Sample_Group%in%c("T1","T2")],levels=c("T2","T1"))
library(ChAMP)
library(stringr)
set.seed(123)
DMR_01_T1_to_T2 <- champ.DMR(beta=as.matrix(myCombat_01[,SampleSheet_01$Sample_Group%in%c("T1","T2")]),
                             pheno=pheno_01_T1_to_T2,
                             method="Bumphunter",
                             compare.group = NULL,
                             arraytype = "EPIC")
dim(DMR_01_T1_to_T2$BumphunterDMR)
write.csv(DMR_01_T1_to_T2$BumphunterDMR, "./DMR_01_T1_to_T2.csv")
index <- apply(DMR_01_T1_to_T2$BumphunterDMR, 1, function(x) which(df_DMP_01_T1_to_T2$CHR ==
                                                gsub("chr","",x[1]) & df_DMP_01_T1_to_T2$MAPINFO >= as.numeric(x[2]) &
                                                df_DMP_01_T1_to_T2$MAPINFO <= as.numeric(x[3])))
DMR_CpG_list_01 <- lapply(index, function(x) df_DMP_01_T1_to_T2[x,])
index1 <- which(mapply(nrow,DMR_CpG_list_01)>=3)
DMR_CpG_table_01 <- do.call(rbind, DMR_CpG_list_01[index1])
dim(DMR_CpG_table_01)
DMR_CpG_table_01 $ DMRindex <- str_split_fixed(rownames(DMR_CpG_table_01), "[.]", 2)[,1]
rownames(DMR_CpG_table_01) <- str_split_fixed(rownames(DMR_CpG_table_01), "[.]", 2)[,2]
saveRDS(DMR_CpG_table_01, "./DMR_CpG_table_01_T1_to_T2.Rds")

pheno_02_T1_to_T2 <- factor(SampleSheet_02$Sample_Group[SampleSheet_02$Sample_Group%in%c("T1","T2")],levels=c("T2","T1"))
set.seed(123)
DMR_02_T1_to_T2 <- champ.DMR(beta=as.matrix(myCombat_02[,SampleSheet_02$Sample_Group%in%c("T1","T2")]),
                             pheno=pheno_02_T1_to_T2,
                             method="Bumphunter",
                             arraytype = "EPIC")

write.csv(DMR_02_T1_to_T2$BumphunterDMR, "./DMR_02_T1_to_T2.csv")
index <- apply(DMR_02_T1_to_T2$BumphunterDMR, 1, function(x) which(df_DMP_02_T1_to_T2$CHR ==
                                                                     gsub("chr","",x[1]) & df_DMP_02_T1_to_T2$MAPINFO >= as.numeric(x[2]) &
                                                                     df_DMP_02_T1_to_T2$MAPINFO <= as.numeric(x[3])))
DMR_CpG_list_02 <- lapply(index, function(x) df_DMP_02_T1_to_T2[x,])
index1 <- which(mapply(nrow,DMR_CpG_list_02)>=3)
DMR_CpG_table_02 <- do.call(rbind, DMR_CpG_list_02[index1])
dim(DMR_CpG_table_02)
DMR_CpG_table_02 $ DMRindex <- str_split_fixed(rownames(DMR_CpG_table_02), "[.]", 2)[,1]
rownames(DMR_CpG_table_02) <- str_split_fixed(rownames(DMR_CpG_table_02), "[.]", 2)[,2]
saveRDS(DMR_CpG_table_02, "./DMR_CpG_table_02_T1_to_T2.Rds")


contained_DMR1 <- data.frame()
contained_DMR2 <- data.frame()

# 检查每对 DMR 是否存在包含关系
DMR_data1 <- DMR_01_T1_to_T2$BumphunterDMR
DMR_data2 <- DMR_02_T1_to_T2$BumphunterDMR
for (i in 1:nrow(DMR_01_T1_to_T2$BumphunterDMR)) {
  print(i)
  for (j in 1:nrow(DMR_02_T1_to_T2$BumphunterDMR)) {
    if (i != j & DMR_data1$seqnames[i] == DMR_data2$seqnames[j]) {
      if (DMR_data1$start[i] >= DMR_data2$start[j] && DMR_data1$end[i] <= DMR_data2$end[j]) {
        contained_DMR2 <- rbind(contained_DMR2, DMR_data2[j, ])
      }
      if (DMR_data2$start[j] >= DMR_data1$start[i] && DMR_data2$end[j] <= DMR_data1$end[i]) {
        contained_DMR1 <- rbind(contained_DMR1, DMR_data1[i, ])
      }
    }
  }
}

contained_DMR1$position <- paste(contained_DMR1$seqnames, contained_DMR1$start, contained_DMR1$end, sep = "_")
contained_DMR2$position <- paste(contained_DMR2$seqnames, contained_DMR2$start, contained_DMR2$end, sep = "_")
dim(contained_DMR2)
contained_DMR2_CpG_table <- DMR_CpG_table_02[DMR_CpG_table_02$DMRindex%in%rownames(contained_DMR2), ]
contained_DMR2_gene_down <- contained_DMR2_CpG_table$gene[contained_DMR2_CpG_table$change_deltaBeta0.05=="down"]
contained_DMR2_gene_down <- contained_DMR2_gene_down[contained_DMR2_gene_down!=""]


go_enrich_res <- go_enrich(as.character(contained_DMR2_gene_down), gene_type="contained_DMR02_gene_down", outdir="./enrich")
kegg_enrich_res <- kegg_enrich(as.character(contained_DMR2_gene_down), gene_type="contained_DMR02_gene_down", outdir="./enrich")
write.table(rbind(go_enrich_res,kegg_enrich_res),file="./enrich/contained_DMR_180-1vs90_gene_down.txt", sep="\t", quote=F)

####################################################################
pheno_03_T1_to_T2 <- factor(SampleSheet_03$Sample_Group[SampleSheet_03$Sample_Group%in%c("T1","T2")],levels=c("T2","T1"))
set.seed(123)
DMR_03_T1_to_T2 <- champ.DMR(beta=as.matrix(myCombat_03[,SampleSheet_03$Sample_Group%in%c("T1","T2")]),
                             pheno=pheno_03_T1_to_T2,
                             method="Bumphunter",
                             arraytype = "EPIC")

write.csv(DMR_03_T1_to_T2$BumphunterDMR, "./DMR_03_T1_to_T2.csv")
index <- apply(DMR_03_T1_to_T2$BumphunterDMR, 1, function(x) which(df_DMP_03_T1_to_T2$CHR ==
                                                                     gsub("chr","",x[1]) & df_DMP_03_T1_to_T2$MAPINFO >= as.numeric(x[2]) &
                                                                     df_DMP_03_T1_to_T2$MAPINFO <= as.numeric(x[3])))
DMR_CpG_list_03 <- lapply(index, function(x) df_DMP_03_T1_to_T2[x,])
index1 <- which(mapply(nrow,DMR_CpG_list_03)>=3)
DMR_CpG_table_03 <- do.call(rbind, DMR_CpG_list_03[index1])
dim(DMR_CpG_table_03)
DMR_CpG_table_03 $ DMRindex <- str_split_fixed(rownames(DMR_CpG_table_03), "[.]", 2)[,1]
rownames(DMR_CpG_table_03) <- str_split_fixed(rownames(DMR_CpG_table_03), "[.]", 2)[,2]
saveRDS(DMR_CpG_table_03, "./DMR_CpG_table_03_T1_to_T2.Rds")

contained_DMR1 <- data.frame()
contained_DMR3 <- data.frame()

DMR_data1 <- read.csv("./DMR_01_T1_to_T2.csv", row.names = 1)
DMR_data2 <- read.csv("./DMR_03_T1_to_T2.csv", row.names = 1)

for (i in 1:nrow(DMR_data1)) {
  print(i)
  for (j in 1:nrow(DMR_data2)) {
    if (i != j & DMR_data1$seqnames[i] == DMR_data2$seqnames[j]) {
      if (DMR_data1$start[i] >= DMR_data2$start[j] && DMR_data1$end[i] <= DMR_data2$end[j]) {
        contained_DMR3 <- rbind(contained_DMR3, DMR_data2[j, ])
      }
      if (DMR_data2$start[j] >= DMR_data1$start[i] && DMR_data2$end[j] <= DMR_data1$end[i]) {
        contained_DMR1 <- rbind(contained_DMR1, DMR_data1[i, ])
      }
    }
  }
}

contained_DMR1$position <- paste(contained_DMR1$seqnames, contained_DMR1$start, contained_DMR1$end, sep = "_")
contained_DMR3$position <- paste(contained_DMR3$seqnames, contained_DMR3$start, contained_DMR3$end, sep = "_")
dim(contained_DMR3)
dim(contained_DMR1)

DMR_CpG_table_03 <- readRDS("./DMR_CpG_table_03_T1_to_T2.Rds")

contained_DMR3_CpG_table <- DMR_CpG_table_03[DMR_CpG_table_03$DMRindex%in%rownames(contained_DMR3), ]
contained_DMR3_gene_down <- contained_DMR3_CpG_table$gene[contained_DMR3_CpG_table$change_deltaBeta0.05=="down"]
contained_DMR3_gene_down <- unique(contained_DMR3_gene_down[contained_DMR3_gene_down!=""])
length(contained_DMR3_gene_down)


go_enrich_res <- go_enrich(as.character(contained_DMR3_gene_down), gene_type="contained_DMR03_gene_down", outdir="./enrich")
kegg_enrich_res <- kegg_enrich(as.character(contained_DMR3_gene_down), gene_type="contained_DMR03_gene_down", outdir="./enrich")
write.table(rbind(go_enrich_res,kegg_enrich_res),file="./enrich/contained_DMR_180-2vs90_gene_down_T2vsT1.txt", sep="\t", quote=F)
