
mydir <- "D:/PhD/major_program/analysis/meth/stp2_DMP/DMR/data"
setwd(mydir)
    
experiment_groups <- rbind(c("01","02"),c("01","03"),c("03","02"))
time_groups <- rbind(c("T1","T2"),c("T2","T3"),c("T1","T3"))
    
library(dplyr)    
res_df <- c()

for (e in 1:nrow(experiment_groups)) {
  experiment_group <- experiment_groups[e,]
  print(experiment_group)
  for (t in 1:nrow(time_groups)) {
      time_group <- time_groups[t,]
      print(time_group)
      contained_DMR1 <- data.frame()
      contained_DMR2 <- data.frame()
      Be_contained_DMR1 <- data.frame()
      Be_contained_DMR2 <- data.frame()
      
      # Check each pair of DMRs for inclusionary relationships
      DMR_data1 <- read.csv(paste0("./DMR_",experiment_group[1],"_",time_group[1],"_to_",time_group[2],".csv"), row.names = 1)
      DMR_data2 <- read.csv(paste0("./DMR_",experiment_group[2],"_",time_group[1],"_to_",time_group[2],".csv"), row.names = 1)
      print(table(DMR_data1$value>0))
      print(table(DMR_data2$value>0))
      
      for (i in 1:nrow(DMR_data1)) {
        for (j in 1:nrow(DMR_data2)) {
          if (i != j & DMR_data1$seqnames[i] == DMR_data2$seqnames[j]) {
            if (DMR_data1$start[i] >= DMR_data2$start[j] && DMR_data1$end[i] <= DMR_data2$end[j]) {
              contained_DMR2 <- rbind(contained_DMR2, DMR_data2[j, ])
              Be_contained_DMR1 <- rbind(Be_contained_DMR1, DMR_data1[i, ])
            }
            if (DMR_data2$start[j] >= DMR_data1$start[i] && DMR_data2$end[j] <= DMR_data1$end[i]) {
              contained_DMR1 <- rbind(contained_DMR1, DMR_data1[i, ])
              Be_contained_DMR2 <- rbind(Be_contained_DMR2, DMR_data2[j, ])
            }
          }
        }
      }
      
      contained_DMR1$position <- paste(contained_DMR1$seqnames, contained_DMR1$start, contained_DMR1$end, sep = "_")
      contained_DMR2$position <- paste(contained_DMR2$seqnames, contained_DMR2$start, contained_DMR2$end, sep = "_")
      Be_contained_DMR1$position <- paste(Be_contained_DMR1$seqnames, Be_contained_DMR1$start, Be_contained_DMR1$end, sep = "_")
      Be_contained_DMR2$position <- paste(Be_contained_DMR2$seqnames, Be_contained_DMR2$start, Be_contained_DMR2$end, sep = "_")
      
      contained_DMR <- as.data.frame(rbind(contained_DMR1,contained_DMR2))
      dim(contained_DMR)
      rownames(contained_DMR) <- paste(c(rep(experiment_group[1],nrow(contained_DMR1)), rep(experiment_group[2],nrow(contained_DMR2))), rownames(contained_DMR), sep="_")
      contained_DMR$group <- c(rep(paste0(experiment_group[1], "contain", experiment_group[2]),nrow(contained_DMR1)), rep(paste0(experiment_group[2], "contain", experiment_group[1]),nrow(contained_DMR2)))
      contained_DMR$DMR_ID <- rownames(contained_DMR)
      
      #对position值一样的行（最多两行）进行：value取均值；行名合并；group的值改为overlap
      result_df <- split(contained_DMR, contained_DMR$position)
      result_df <- lapply(result_df, function(df){df <- df[,c(1:6,18:19)];
      if(nrow(df)==2){
        df$group <- "overlap"; df$value <- df$value[which(abs(df$value)==max(abs(df$value)))]; df$DMR_ID <- paste0(df$DMR_ID, collapse = "and")
      }
      return(df)})  
      contained_DMR_integrated <- as.data.frame(unique(do.call(rbind,result_df)))
      rownames(contained_DMR_integrated) <- contained_DMR_integrated$DMR_ID      
      contained_DMR_integrated <- contained_DMR_integrated[,-ncol(contained_DMR_integrated)]
      
      if(nrow(contained_DMR)!=0){
        write.table(contained_DMR,file=paste0("./contained_DMR_", experiment_group[2],"vs",experiment_group[1],"_", time_group[2],"vs",time_group[1],".txt"), sep="\t", quote=F)
        write.table(contained_DMR_integrated,file=paste0("./contained_DMR_integrated_", experiment_group[2],"vs",experiment_group[1],"_", time_group[2],"vs",time_group[1],".txt"), sep="\t", quote=F)
      }
      
      DMR_data1$position <- paste(DMR_data1$seqnames, DMR_data1$start, DMR_data1$end, sep = "_")
      DMR_data2$position <- paste(DMR_data2$seqnames, DMR_data2$start, DMR_data2$end, sep = "_")
      DMR_data1_unique <-  DMR_data1[!DMR_data1$position%in%c(contained_DMR1$position,Be_contained_DMR1$position),]
      DMR_data2_unique <-  DMR_data2[!DMR_data2$position%in%c(contained_DMR2$position,Be_contained_DMR2$position),]
      write.table(DMR_data1_unique,file=paste0("./unique_DMR_for",experiment_group[1],"_", experiment_group[2],"vs",experiment_group[1],"_", time_group[2],"vs",time_group[1],".txt"), sep="\t", quote=F)
      write.table(DMR_data2_unique,file=paste0("./unique_DMR_for",experiment_group[2],"_", experiment_group[2],"vs",experiment_group[1],"_", time_group[2],"vs",time_group[1],".txt"), sep="\t", quote=F)
      
      
      DMR_CpG_table2 <- readRDS(paste0(".DMR_CpG_table_",experiment_group[2],"_",time_group[1],"_to_",time_group[2],".Rds"))
      DMR_CpG_table1 <- readRDS(paste0(".DMR_CpG_table_",experiment_group[1],"_",time_group[1],"_to_",time_group[2],".Rds"))
      
      contained_DMR_CpG_table <- rbind(DMR_CpG_table2[DMR_CpG_table2$DMRindex%in%rownames(contained_DMR2), ], DMR_CpG_table1[DMR_CpG_table1$DMRindex%in%rownames(contained_DMR1), ])
      contained_DMR_gene_down <- contained_DMR_CpG_table$gene[contained_DMR_CpG_table$change_deltaBeta0.05=="down"]
      contained_DMR_gene_down <- unique(contained_DMR_gene_down[contained_DMR_gene_down!=""])
      contained_DMR_gene_up <- contained_DMR_CpG_table$gene[contained_DMR_CpG_table$change_deltaBeta0.05=="up"]
      contained_DMR_gene_up <- unique(contained_DMR_gene_up[contained_DMR_gene_up!=""])
      
      go_enrich_res_down <- c()
      kegg_enrich_res_down <- c()
      if(length(contained_DMR_gene_down)>1){
        go_enrich_res_down <- go_enrich(as.character(contained_DMR_gene_down), gene_type="contained_DMR02_gene_down", outdir="./enrich")
        kegg_enrich_res_down <- kegg_enrich(as.character(contained_DMR_gene_down), gene_type="contained_DMR02_gene_down", outdir="./enrich")
        if(!is.null(go_enrich_res_down)|!is.null(kegg_enrich_res_down)){
          write.table(rbind(go_enrich_res_down,kegg_enrich_res_down),file=paste0("./enrich/contained_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_down.txt"), sep="\t", quote=F)
        }
      }
      
      go_enrich_res_up <- c()
      kegg_enrich_res_up <- c()
      if(length(contained_DMR_gene_up)>1){
        go_enrich_res_up <- go_enrich(as.character(contained_DMR_gene_up), gene_type="contained_DMR02_gene_up", outdir="./enrich")
        kegg_enrich_res_up <- kegg_enrich(as.character(contained_DMR_gene_up), gene_type="contained_DMR02_gene_up", outdir="./enrich")
        if(!is.null(go_enrich_res_up)|!is.null(kegg_enrich_res_up)){
          write.table(rbind(go_enrich_res_up,kegg_enrich_res_up),file=paste0("./enrich/contained_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_up.txt"), sep="\t", quote=F)
        }
      }
      
      
      #Enrichment for DMR_data1_unique
      DMR_data1_unique_CpG_table <- DMR_CpG_table1[DMR_CpG_table1$DMRindex%in%rownames(DMR_data1_unique), ]
      DMR_data1_unique_gene_down <- DMR_data1_unique_CpG_table$gene[DMR_data1_unique_CpG_table$change_deltaBeta0.05=="down"]
      DMR_data1_unique_gene_down <- unique(DMR_data1_unique_gene_down[DMR_data1_unique_gene_down!=""])
      DMR_data1_unique_gene_up <- DMR_data1_unique_CpG_table$gene[DMR_data1_unique_CpG_table$change_deltaBeta0.05=="up"]
      DMR_data1_unique_gene_up <- unique(DMR_data1_unique_gene_up[DMR_data1_unique_gene_up!=""])
      
      DMR_data1_go_enrich_res_down <- c()
      DMR_data1_kegg_enrich_res_down <- c()
      if(length(DMR_data1_unique_gene_down)>1){
        DMR_data1_go_enrich_res_down <- go_enrich(as.character(DMR_data1_unique_gene_down), gene_type="DMR_data1_unique02_gene_down", outdir="./enrich")
        DMR_data1_kegg_enrich_res_down <- kegg_enrich(as.character(DMR_data1_unique_gene_down), gene_type="DMR_data1_unique02_gene_down", outdir="./enrich")
        if(!is.null(DMR_data1_go_enrich_res_down)|!is.null(DMR_data1_kegg_enrich_res_down)){
          write.table(rbind(DMR_data1_go_enrich_res_down,DMR_data1_kegg_enrich_res_down),file=paste0("./enrich/",experiment_group[1],"_unique_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_down.txt"), sep="\t", quote=F)
        }
      }
      
      DMR_data1_go_enrich_res_up <- c()
      DMR_data1_kegg_enrich_res_up <- c()
      if(length(DMR_data1_unique_gene_up)>1){
        DMR_data1_go_enrich_res_up <- go_enrich(as.character(DMR_data1_unique_gene_up), gene_type="DMR_data1_unique02_gene_up", outdir="./enrich")
        DMR_data1_kegg_enrich_res_up <- kegg_enrich(as.character(DMR_data1_unique_gene_up), gene_type="DMR_data1_unique02_gene_up", outdir="./enrich")
        if(!is.null(DMR_data1_go_enrich_res_up)|!is.null(DMR_data1_kegg_enrich_res_up)){
          write.table(rbind(DMR_data1_go_enrich_res_up,DMR_data1_kegg_enrich_res_up),file=paste0("./enrich/",experiment_group[1],"_unique_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_up.txt"), sep="\t", quote=F)
        }
      }
      
      #Enrichment for DMR_data2_unique
      DMR_data2_unique_CpG_table <- DMR_CpG_table2[DMR_CpG_table2$DMRindex%in%rownames(DMR_data2_unique), ]
      DMR_data2_unique_gene_down <- DMR_data2_unique_CpG_table$gene[DMR_data2_unique_CpG_table$change_deltaBeta0.05=="down"]
      DMR_data2_unique_gene_down <- unique(DMR_data2_unique_gene_down[DMR_data2_unique_gene_down!=""])
      DMR_data2_unique_gene_up <- DMR_data2_unique_CpG_table$gene[DMR_data2_unique_CpG_table$change_deltaBeta0.05=="up"]
      DMR_data2_unique_gene_up <- unique(DMR_data2_unique_gene_up[DMR_data2_unique_gene_up!=""])
      
      DMR_data2_go_enrich_res_down <- c()
      DMR_data2_kegg_enrich_res_down <- c()
      if(length(DMR_data2_unique_gene_down)>1){
        DMR_data2_go_enrich_res_down <- go_enrich(as.character(DMR_data2_unique_gene_down), gene_type="DMR_data2_unique02_gene_down", outdir="./enrich")
        DMR_data2_kegg_enrich_res_down <- kegg_enrich(as.character(DMR_data2_unique_gene_down), gene_type="DMR_data2_unique02_gene_down", outdir="./enrich")
        if(!is.null(DMR_data2_go_enrich_res_down)|!is.null(DMR_data2_kegg_enrich_res_down)){
          write.table(rbind(DMR_data2_go_enrich_res_down,DMR_data2_kegg_enrich_res_down),file=paste0("./enrich/",experiment_group[2],"_unique_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_down.txt"), sep="\t", quote=F)
        }
      }
      
      DMR_data2_go_enrich_res_up <- c()
      DMR_data2_kegg_enrich_res_up <- c()
      if(length(DMR_data2_unique_gene_up)>1){
        DMR_data2_go_enrich_res_up <- go_enrich(as.character(DMR_data2_unique_gene_up), gene_type="DMR_data2_unique02_gene_up", outdir="./enrich")
        DMR_data2_kegg_enrich_res_up <- kegg_enrich(as.character(DMR_data2_unique_gene_up), gene_type="DMR_data2_unique02_gene_up", outdir="./enrich")
        if(!is.null(DMR_data2_go_enrich_res_up)|!is.null(DMR_data2_kegg_enrich_res_up)){
          write.table(rbind(DMR_data2_go_enrich_res_up,DMR_data2_kegg_enrich_res_up),file=paste0("./enrich/",experiment_group[2],"_unique_DMR_", experiment_group[2],"vs",experiment_group[1],"_gene_up.txt"), sep="\t", quote=F)
        }
      }
      DMR_summary_fun <- function(x){
        tab <- table(x$value>=0)
        if(is.na(tab["TRUE"])){tab["TRUE"] <- 0}
        if(is.na(tab["FALSE"])){tab["FALSE"] <- 0}
        return(c(tab["TRUE"],tab["FALSE"]))
      } 

      res <- c(paste0(experiment_group[2], "vs", experiment_group[1]), paste0(time_group[2], "vs", time_group[1]), DMR_summary_fun(DMR_data2), DMR_summary_fun(DMR_data1), DMR_summary_fun(contained_DMR1), DMR_summary_fun(contained_DMR2), length(intersect(contained_DMR1$position, contained_DMR2$position)),
                        DMR_summary_fun(DMR_data2_unique), DMR_summary_fun(DMR_data1_unique),
                        length(contained_DMR_gene_down), length(contained_DMR_gene_up), down_gene_fun=paste0(go_enrich_res_down$Description, collapse=","), up_gene_fun=paste0(go_enrich_res_up$Description, collapse=","),
                        length(DMR_data1_unique_gene_down), length(DMR_data1_unique_gene_up), down_gene_fun=paste0(DMR_data1_go_enrich_res_down$Description, collapse=","), up_gene_fun=paste0(DMR_data1_go_enrich_res_up$Description, collapse=","),
                        length(DMR_data2_unique_gene_down), length(DMR_data2_unique_gene_up), down_gene_fun=paste0(DMR_data2_go_enrich_res_down$Description, collapse=","), up_gene_fun=paste0(DMR_data2_go_enrich_res_up$Description, collapse=","))
      res_df <- rbind(res_df, res)
  }
}

write.csv(res_df, ".DMR_intersection_new.csv", row.names = F)
