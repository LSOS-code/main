##### description: Statistics on the number of DMPs and presentation of the changing trend of DMPs in T2 vs T1

#load packages
library(RColorBrewer)
library(ggplot2)

##### 01Start: Statistics on the number of hypermethylated and hypomethylated DMPs #####
df_DMP_list <- readRDS("Differential methylation analysis/result/DMP_res_M90_M180-1.Rds")
df_DMP_01_T1_to_T2 <- df_DMP_list[[1]]
df_DMP_01_T2_to_T3 <- df_DMP_list[[2]]
df_DMP_01_T1_to_T3 <- df_DMP_list[[3]]

df_DMP_02_T1_to_T2 <- df_DMP_list[[4]]
df_DMP_02_T2_to_T3 <- df_DMP_list[[5]]
df_DMP_02_T1_to_T3 <- df_DMP_list[[6]]

DMP_list_01_deltaBeta0.05 <- list(`up_T1_to_T2` = rownames(df_DMP_01_T1_to_T2[df_DMP_01_T1_to_T2$change_deltaBeta0.05=="up",]),
                    `down_T1_to_T2` = rownames(df_DMP_01_T1_to_T2[df_DMP_01_T1_to_T2$change_deltaBeta0.05=="down",]),
					`stable_T1_to_T2` = rownames(df_DMP_01_T1_to_T2[df_DMP_01_T1_to_T2$change_deltaBeta0.05=="stable",]),
                    `up_T2_to_T3` = rownames(df_DMP_01_T2_to_T3[df_DMP_01_T2_to_T3$change_deltaBeta0.05=="up",]),
                    `down_T2_to_T3` = rownames(df_DMP_01_T2_to_T3[df_DMP_01_T2_to_T3$change_deltaBeta0.05=="down",]),
					`stable_T2_to_T3` = rownames(df_DMP_01_T2_to_T3[df_DMP_01_T2_to_T3$change_deltaBeta0.05=="stable",]))
lapply(DMP_list_01_deltaBeta0.05, length)

DMP_list_02_deltaBeta0.05 <- list(`up_T1_to_T2` = rownames(df_DMP_02_T1_to_T2[df_DMP_02_T1_to_T2$change_deltaBeta0.05=="up",]),
                    `down_T1_to_T2` = rownames(df_DMP_02_T1_to_T2[df_DMP_02_T1_to_T2$change_deltaBeta0.05=="down",]),
					`stable_T1_to_T2` = rownames(df_DMP_02_T1_to_T2[df_DMP_02_T1_to_T2$change_deltaBeta0.05=="stable",]),
                    `up_T2_to_T3` = rownames(df_DMP_02_T2_to_T3[df_DMP_02_T2_to_T3$change_deltaBeta0.05=="up",]),
                    `down_T2_to_T3` = rownames(df_DMP_02_T2_to_T3[df_DMP_02_T2_to_T3$change_deltaBeta0.05=="down",]),
					`stable_T2_to_T3` = rownames(df_DMP_02_T2_to_T3[df_DMP_02_T2_to_T3$change_deltaBeta0.05=="stable",]))
lapply(DMP_list_02_deltaBeta0.05, length)

DMP_num_sta_0.05 <- mapply(function(x) table(x$change_deltaBeta0.05), df_DMP_list)
colnames(DMP_num_sta_0.05) <- c("01_T1_to_T2", "01_T2_to_T3", "02_T1_to_T2", "02_T2_to_T3")
##### 01End: Statistics on the number of hypermethylated and hypomethylated DMPs #####

##### 02Start: Count the number of DMPs and categorise them according to the change process in the T2_vs_T1 and T3_vs_T2 #####
# For M90
state <- c("up", "stable", "down")
state_groups <- unique(t(combn(c(1:3,1:3),2)))
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")


DMP_list_trends <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_genes <- intersect(DMP_list_01_deltaBeta0.05[[group1]], DMP_list_01_deltaBeta0.05[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(DMP_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)
mapply(length, DMP_list_trends )
DMP_list_trends <- DMP_list_trends[-grep("stable_stable", names(DMP_list_trends))]
DMP_list_trends <- DMP_list_trends[mapply(length, DMP_list_trends )!=0]

saveRDS(DMP_list_trends, file = "Differential methylation analysis/result/90-day spaceflight_M90/DMP_list_trend.Rds")

DMP_list_trend_01 <- readRDS("Differential methylation analysis/result/90-day spaceflight_M90/DMP_list_trend.Rds")
DMP_num_01 <- mapply(length, DMP_list_trend_01)
DMP_num_df_01 <- data.frame(trend=names(DMP_num_01), number=DMP_num_01)
DMP_num_df_01$T1_to_T2 <- stringr::str_split_i(DMP_num_df_01$trend, "_", 1)
DMP_num_df_01$T2_to_T3 <- stringr::str_split_i(DMP_num_df_01$trend, "_", 2)
write.table(DMP_num_df_01, file = "Differential methylation analysis/result/90-day spaceflight_M90/DMP_number_statistic.txt", sep = "\t", quote = F, row.names = F)

DMP_num_df_01 <- DMP_num_df_01[DMP_num_df_01$T1_to_T2!="stable",]
DMP_num_df_01$T1_to_T2<-factor(DMP_num_df_01$T1_to_T2,
                       levels = c("up","stable","down"))
DMP_num_df_01$T2_to_T3<-factor(DMP_num_df_01$T2_to_T3,
                               levels = c("up","stable","down"))
#Creating an alluvium map	
ggplot(data=DMP_num_df_01,
       aes(axis1=T1_to_T2,axis2=T2_to_T3,
           y=number))+
  geom_alluvium(aes(fill=T1_to_T2),
                #size=3,
                #color="white",
                width = 0.1,
                aes.bind = "flows")+
  geom_stratum(fill=c("#73B9CE", "#C1548E","#73B9CE", "darkgray", "#C1548E"),
               #color="white",
               #size=3,
               width=0.1)+
  scale_fill_manual(breaks = c("up","stable","down"),
                    values = c("#C1548E", "darkgray", "#73B9CE"),
                    labels = c("up","stable","down")) +
  scale_color_manual(breaks = c("up","stable","down"),
                     values = c("#C1548E", "darkgray", "#73B9CE"),
                     labels = c("up","stable","down")) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_continuous(breaks = c(1,2),
                     labels = c("T1_to_T2","T2_to_T3"),
                     expand = expansion(mult = c(0,0)))+
  coord_cartesian(clip="off") +
  theme_classic()
  
#For M180-1
state <- c("up", "stable", "down")
state_groups <- unique(t(combn(c(1:3,1:3),2)))
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")


DMP_list_trends <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_genes <- intersect(DMP_list_02_deltaBeta0.05[[group1]], DMP_list_02_deltaBeta0.05[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(DMP_list_trends) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)
mapply(length, DMP_list_trends )
DMP_list_trends <- DMP_list_trends[-grep("stable_stable", names(DMP_list_trends))]
DMP_list_trends <- DMP_list_trends[mapply(length, DMP_list_trends )!=0]

saveRDS(DMP_list_trends, file = "Differential methylation analysis/result/180-day spaceflight_M180-1/DMP_list_trend.Rds")

DMP_list_trend_02 <- readRDS("Differential methylation analysis/result/180-day spaceflight_M180-1/DMP_list_trend.Rds")
DMP_num_02 <- mapply(length, DMP_list_trend_02)
DMP_num_df_02 <- data.frame(trend=names(DMP_num_02), number=DMP_num_02)
DMP_num_df_02$T1_to_T2 <- stringr::str_split_i(DMP_num_df_02$trend, "_", 1)
DMP_num_df_02$T2_to_T3 <- stringr::str_split_i(DMP_num_df_02$trend, "_", 2)
write.table(DMP_num_df_02, file = "Differential methylation analysis/result/180-day spaceflight_M180-1/DMP_number_statistic.txt", sep = "\t", quote = F, row.names = F)

DMP_num_df_02 <- DMP_num_df_02[DMP_num_df_02$T1_to_T2!="stable",]
DMP_num_df_02$T1_to_T2<-factor(DMP_num_df_02$T1_to_T2,
                       levels = c("up","stable","down"))
DMP_num_df_02$T2_to_T3<-factor(DMP_num_df_02$T2_to_T3,
                               levels = c("up","stable","down"))
#Creating an alluvium map	
ggplot(data=DMP_num_df_02,
       aes(axis1=T1_to_T2,axis2=T2_to_T3,
           y=number))+
  geom_alluvium(aes(fill=T1_to_T2),
                #size=3,
                #color="white",
                width = 0.1,
                aes.bind = "flows")+
  geom_stratum(fill=c("#73B9CE", "#C1548E","#73B9CE", "darkgray", "#C1548E"),
               #color="white",
               #size=3,
               width=0.1)+
  scale_fill_manual(breaks = c("up","stable","down"),
                    values = c("#C1548E", "darkgray", "#73B9CE"),
                    labels = c("up","stable","down")) +
  scale_color_manual(breaks = c("up","stable","down"),
                     values = c("#C1548E", "darkgray", "#73B9CE"),
                     labels = c("up","stable","down")) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))+
  scale_x_continuous(breaks = c(1,2),
                     labels = c("T1_to_T2","T2_to_T3"),
                     expand = expansion(mult = c(0,0)))+
  coord_cartesian(clip="off") +
  theme_classic()
##### 02End: Count the number of DMPs and categorise them according to the change process in the T2_vs_T1 and T3_vs_T2 #####


##### 03Start: Trend plot for differentially methylated probes #####

#trend_matrix: col:gene; rownames:time by order; value:mean of methlation beta value of the genes.
trend_plot <- function(trend_matrix, my_col, title){
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
    pdata<-data.frame(x,lower = quantile_value[k,],upper = quantile_value[101-k,])
    pdata.list[[k]]<-pdata
  }
  myPalette <- colorRampPalette(my_col)(45)
  plot.trend<-ggplot()+
    geom_rect(aes(xmin = -Inf,xmax = 1,ymin = -Inf,ymax = Inf),fill = "white",alpha = 0.4)+
    geom_rect(aes(xmin = 1,xmax = length(time),ymin = -Inf,ymax = Inf),fill = "#DEF6F3",alpha = 0.4)+
    geom_rect(aes(xmin = length(time),xmax = Inf,ymin = -Inf,ymax = Inf),fill = "white",alpha = 0.4)+
    theme(panel.grid=element_blank(),panel.border=element_rect(fill=NA,color="black", size=1))
  
  for(k in 6:50){
    pdata<-pdata.list[[k]]
    plot.trend<-plot.trend+geom_ribbon(data = pdata,aes(ymin=lower, ymax=upper, x=x), fill = myPalette[k-5], alpha = 1)
  }
  plot.trend<-plot.trend+
    theme_bw()+
    xlab("Time")+
    ylab("deltaBeta(vs T1)")+
    labs(title = title)+
    theme(panel.grid =element_blank()) + 
    scale_x_continuous(breaks=seq(1, length(time), 1), labels = c("T1","T2","T3"))
  return(plot.trend)
}

#For M90
probes <- rownames(df_DMP_01_T1_to_T2)
deltaBeta_df_01 <- data.frame(probe_name=probes, T1_value=0, T2_value=df_DMP_01_T1_to_T2$deltaBeta,
                          T3_value=df_DMP_01_T1_to_T3[probes, "deltaBeta"], change_type=df_DMP_01_T1_to_T2$change_deltaBeta0.05)
write.csv(deltaBeta_df_01, file = "Changing trend of DMPs M90.csv", row.names=F)
deltaBeta_df_01_up <- deltaBeta_df_01[deltaBeta_df_01$change_type=="up",2:4]
my_col <- c("#F0D8E5", "#C1548E")
up_plot_1 <- trend_plot(deltaBeta_df_01_up, my_col, title = "Hyper-methylated probes in T2 vs T1")
deltaBeta_df_01_down <- deltaBeta_df_01[deltaBeta_df_01$change_type=="down",2:4]
my_col <- c("#C6E4F5", "#6BACD1")
down_plot_1 <- trend_plot(deltaBeta_df_01_down, my_col, title = "Hypo-methylated probes in T2 vs T1")

##For M180-1
probes <- rownames(df_DMP_01_T1_to_T2)
deltaBeta_df_02 <- data.frame(probe_name=probes, T1_value=0, T2_value=df_DMP_02_T1_to_T2$deltaBeta,
                          T3_value=df_DMP_02_T1_to_T3[probes, "deltaBeta"], change_type=df_DMP_02_T1_to_T2$change_deltaBeta0.05)
write.csv(deltaBeta_df_02, file = "Changing trend of DMPs M180-1.csv", row.names=F)
deltaBeta_df_02_up <- deltaBeta_df_02[deltaBeta_df_02$change_type=="up",2:4]
my_col <- c("#F0D8E5", "#C1548E")
up_plot_2 <- trend_plot(deltaBeta_df_02_up, my_col, title = "Hyper-methylated probes in T2 vs T1")
deltaBeta_df_02_down <- deltaBeta_df_02[deltaBeta_df_02$change_type=="down",2:4]
my_col <- c("#C6E4F5", "#6BACD1")
down_plot_2 <- trend_plot(deltaBeta_df_02_down, my_col, title = "Hypo-methylated probes in T2 vs T1")
ggpubr::ggarrange(up_plot_1, down_plot_1, up_plot_2, down_plot_2, ncol = 2, nrow = 2)
##### 03End: Trend plot for differentially methylated probes #####


##### 04Start: Methylation changing trend of genes #####
#Count the overall methylation changes in the promoter region of a gene, only if all probes are consistently hypermethylated or hypomethylated will the gene be labelled up or down accordingly
#For M90
state <- c("up", "stable", "down")
state_groups <- unique(t(combn(c(1:3,1:3),2)))
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")
trends <- trends[!apply(trends, 1, function(x){x[1]==x[2] & x[2] == "stable"}), ]

DMP_trends_list_01 <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_probes <- intersect(DMP_list_01_deltaBeta0.05[[group1]], DMP_list_01_deltaBeta0.05[[group2]])
  print(length(trend_probes))
  return(trend_probes)
})
names(DMP_trends_list_01) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)

DMG_list_01_deltaBeta0.05 <- DMP_list_01_deltaBeta0.05
for(i in 1:length(DMP_list_01_deltaBeta0.05)){
  probe_CpG_info<-df_DMP_01_T1_to_T2[DMP_list_01_deltaBeta0.05[[i]], ]
  promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")
  DM_genes<-as.character(probe_CpG_info[probe_CpG_info$feature%in%promoter_region,"gene"])
  DM_genes<-DM_genes[DM_genes!=""]
  DM_genes<-unique(na.omit(DM_genes))
  DMG_list_01_deltaBeta0.05[[i]]<-DM_genes
}
DMG_list_01_deltaBeta0.05$stable_T1_to_T2 <- DMG_list_01_deltaBeta0.05$stable_T1_to_T2[!DMG_list_01_deltaBeta0.05$stable_T1_to_T2%in%c(DMG_list_01_deltaBeta0.05$up_T1_to_T2, DMG_list_01_deltaBeta0.05$down_T1_to_T2)]
common_genes <- intersect(DMG_list_01_deltaBeta0.05$up_T1_to_T2, DMG_list_01_deltaBeta0.05$down_T1_to_T2)
DMG_list_01_deltaBeta0.05$up_T1_to_T2 <- DMG_list_01_deltaBeta0.05$up_T1_to_T2[!DMG_list_01_deltaBeta0.05$up_T1_to_T2%in%common_genes]
DMG_list_01_deltaBeta0.05$down_T1_to_T2 <- DMG_list_01_deltaBeta0.05$down_T1_to_T2[!DMG_list_01_deltaBeta0.05$down_T1_to_T2%in%common_genes]
DMG_list_01_deltaBeta0.05$stable_T2_to_T3 <- DMG_list_01_deltaBeta0.05$stable_T2_to_T3[!DMG_list_01_deltaBeta0.05$stable_T2_to_T3%in%c(DMG_list_01_deltaBeta0.05$up_T2_to_T3, DMG_list_01_deltaBeta0.05$down_T2_to_T3)]
common_genes <- intersect(DMG_list_01_deltaBeta0.05$up_T2_to_T3, DMG_list_01_deltaBeta0.05$down_T2_to_T3)
DMG_list_01_deltaBeta0.05$up_T2_to_T3 <- DMG_list_01_deltaBeta0.05$up_T2_to_T3[!DMG_list_01_deltaBeta0.05$up_T2_to_T3%in%common_genes]
DMG_list_01_deltaBeta0.05$down_T2_to_T3 <- DMG_list_01_deltaBeta0.05$down_T2_to_T3[!DMG_list_01_deltaBeta0.05$down_T2_to_T3%in%common_genes]

mapply(length, DMG_list_01_deltaBeta0.05)

patterns <- names(DMP_trends_list_01)
pattern_gene_list_01 <- lapply(patterns, function(x){
  x <- unlist(strsplit(x, "_"))
  group1 <- paste0(x[1], "_", "T1_to_T2")
  group2 <- paste0(x[2], "_", "T2_to_T3")
  trend_genes <- intersect(DMG_list_01_deltaBeta0.05[[group1]], DMG_list_01_deltaBeta0.05[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(pattern_gene_list_01) <- patterns
mapply(length, pattern_gene_list_01)
ggvenn::ggvenn(pattern_gene_list_01, c("up_stable", "up_down", "down_up", "down_stable")) 
saveRDS(pattern_gene_list_01, file="Differential methylation analysis/result/180-day spaceflight_M180-1/pattern_gene_list.Rds")

#For M180-1
state <- c("up", "stable", "down")
state_groups <- unique(t(combn(c(1:3,1:3),2)))
trends <- t(apply(state_groups, 1, function(x){return(c(state[x]))}))
trends <- as.data.frame(trends)
colnames(trends) <- c("T1_to_T2","T2_to_T3")
trends <- trends[!apply(trends, 1, function(x){x[1]==x[2] & x[2] == "stable"}), ]

DMP_trends_list_02 <- apply(trends, 1, function(x){
  group1 <- paste0(x[1], "_", colnames(trends)[1])
  group2 <- paste0(x[2], "_", colnames(trends)[2])
  trend_probes <- intersect(DMP_list_02_deltaBeta0.05[[group1]], DMP_list_02_deltaBeta0.05[[group2]])
  print(length(trend_probes))
  return(trend_probes)
})
names(DMP_trends_list_02) <- paste0(trends$T1_to_T2, "_", trends$T2_to_T3)

DMG_list_02_deltaBeta0.05 <- DMP_list_02_deltaBeta0.05
for(i in 1:length(DMP_list_02_deltaBeta0.05)){
  probe_CpG_info<-df_DMP_02_T1_to_T2[DMP_list_02_deltaBeta0.05[[i]], ]
  promoter_region<-c("TSS1500","TSS200","5'UTR","1stExon")
  DM_genes<-as.character(probe_CpG_info[probe_CpG_info$feature%in%promoter_region,"gene"])
  DM_genes<-DM_genes[DM_genes!=""]
  DM_genes<-unique(na.omit(DM_genes))
  DMG_list_02_deltaBeta0.05[[i]]<-DM_genes
}
DMG_list_02_deltaBeta0.05$stable_T1_to_T2 <- DMG_list_02_deltaBeta0.05$stable_T1_to_T2[!DMG_list_02_deltaBeta0.05$stable_T1_to_T2%in%c(DMG_list_02_deltaBeta0.05$up_T1_to_T2, DMG_list_02_deltaBeta0.05$down_T1_to_T2)]
common_genes <- intersect(DMG_list_02_deltaBeta0.05$up_T1_to_T2, DMG_list_02_deltaBeta0.05$down_T1_to_T2)
DMG_list_02_deltaBeta0.05$up_T1_to_T2 <- DMG_list_02_deltaBeta0.05$up_T1_to_T2[!DMG_list_02_deltaBeta0.05$up_T1_to_T2%in%common_genes]
DMG_list_02_deltaBeta0.05$down_T1_to_T2 <- DMG_list_02_deltaBeta0.05$down_T1_to_T2[!DMG_list_02_deltaBeta0.05$down_T1_to_T2%in%common_genes]
DMG_list_02_deltaBeta0.05$stable_T2_to_T3 <- DMG_list_02_deltaBeta0.05$stable_T2_to_T3[!DMG_list_02_deltaBeta0.05$stable_T2_to_T3%in%c(DMG_list_02_deltaBeta0.05$up_T2_to_T3, DMG_list_02_deltaBeta0.05$down_T2_to_T3)]
common_genes <- intersect(DMG_list_02_deltaBeta0.05$up_T2_to_T3, DMG_list_02_deltaBeta0.05$down_T2_to_T3)
DMG_list_02_deltaBeta0.05$up_T2_to_T3 <- DMG_list_02_deltaBeta0.05$up_T2_to_T3[!DMG_list_02_deltaBeta0.05$up_T2_to_T3%in%common_genes]
DMG_list_02_deltaBeta0.05$down_T2_to_T3 <- DMG_list_02_deltaBeta0.05$down_T2_to_T3[!DMG_list_02_deltaBeta0.05$down_T2_to_T3%in%common_genes]

mapply(length, DMG_list_02_deltaBeta0.05)

patterns <- names(DMP_trends_list_02)
pattern_gene_list_02 <- lapply(patterns, function(x){
  x <- unlist(strsplit(x, "_"))
  group1 <- paste0(x[1], "_", "T1_to_T2")
  group2 <- paste0(x[2], "_", "T2_to_T3")
  trend_genes <- intersect(DMG_list_02_deltaBeta0.05[[group1]], DMG_list_02_deltaBeta0.05[[group2]])
  print(length(trend_genes))
  return(trend_genes)
})
names(pattern_gene_list_02) <- patterns
mapply(length, pattern_gene_list_02)
ggvenn::ggvenn(pattern_gene_list_02, c("up_stable", "up_down", "down_up", "down_stable")) 
saveRDS(pattern_gene_list_02, file="Differential methylation analysis/result/180-day spaceflight_M180-1/pattern_gene_list.Rds")
##### 04End: Methylation changing trend of genes #####