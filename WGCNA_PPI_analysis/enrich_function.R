library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library("stringr")

go_enrich <- function(gene_symbols, gene_type){
  if(length(gene_symbols)!=0){
    cols <- c("SYMBOL", "ENSEMBL", 'ENTREZID')
    gene = AnnotationDbi::select(org.Hs.eg.db, keys=gene_symbols, columns=cols, keytype="SYMBOL")
    gene_type <- gene_type
    
    print(length(gene$ENTREZID))
    IDs <- unique(na.omit(gene$ENTREZID))
    ego_BP <- enrichGO(gene = IDs, OrgDb= org.Hs.eg.db, ont = "BP",
                       pAdjustMethod = "BH", minGSSize = 10,
                       pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE)
    ego_BP <- clusterProfiler::simplify(ego_BP)
    go.df_all <- data.frame(ego_BP)
    go.df_all$gene_type <- rep(gene_type, nrow(go.df_all))
    print(dim(go.df_all))
    if(nrow(go.df_all)>0){
      go.df <- go.df_all
      if(nrow(go.df_all)>30){go.df <- go.df_all[1:30,]}
      go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
      go.df$GeneRatio <- unlist(lapply(unlist(go.df$GeneRatio),function(x) eval(parse(text=x))))
      go.df$GeneRatio <- round(go.df$GeneRatio, 3)
      go_bar <- ggplot(data = go.df, 
                       aes(x = Description, y = GeneRatio,fill = pvalue))+ 
        scale_fill_gradient(low = "#E41A1C",high = "white")+
        geom_bar(stat = "identity",width = 0.9)+ 
        coord_flip()+theme_bw()+ 
        scale_x_discrete(labels = function(x) str_wrap(x,width = 80))+ 
        labs(x = "GO terms",y = "GeneRatio",title = paste0("Enriched GO Terms of ", gene_type))+ 
        theme(axis.title = element_text(size = 13), 
              axis.text = element_text(size = 11), 
              #axis.text.y = element_text(color = x.color),
              plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
              legend.title = element_text(size = 13), 
              legend.text = element_text(size = 11), 
              plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
      ggsave(
        filename = paste0("./image/", gene_type, "_Go_barplot.pdf"), 
        plot = go_bar,
        width = 11,             
        height = 7,            
        units = "in",          
        dpi = 300              
      )
      return(go.df_all)
    }
    
  }
}

#R.utils::setOption("clusterProfiler.download.method",'auto') 
kegg_enrich <- function(gene_symbols, gene_type){
  if(length(gene_symbols)!=0){
    cols <- c("ENSEMBL", "SYMBOL", "ENTREZID")
    gene = AnnotationDbi::select(org.Hs.eg.db, keys=gene_symbols, columns=cols, keytype="SYMBOL")
    gene_type <- gene_type
    print(length(gene$ENTREZID))
    
    IDs <- unique(na.omit(gene$ENTREZID))
    kk <- enrichKEGG(gene = IDs, organism ="hsa", pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1, minGSSize = 10,
                     use_internal_data = F)
    kegg.df_all <- data.frame(kk)
    geneSymbols_vec <- lapply(kegg.df_all$geneID, function(x){
      geneIDs <- unlist(strsplit(x, "/"))
      geneSymbols <- gene[match(geneIDs, gene$ENTREZID), "SYMBOL"]
      geneSymbols <- paste0(geneSymbols, collapse = "/")
      return(geneSymbols)
    })
    kegg.df_all$geneID <- unlist(geneSymbols_vec)
    kegg.df_all$gene_type <- rep(gene_type, nrow(kegg.df_all))
    print(dim(kegg.df_all ))
    if(nrow(kegg.df_all)>0){
      kegg.df <- kegg.df_all
      if(nrow(kegg.df_all)>30){kegg.df <- kegg.df_all[1:30,]}
      kegg.df$Description <- factor(kegg.df$Description,levels = rev(kegg.df$Description))
      kegg.df$GeneRatio <- unlist(lapply(unlist(kegg.df$GeneRatio),function(x) eval(parse(text=x))))
      kegg.df$GeneRatio <- round(kegg.df$GeneRatio, 3)
      kegg_bar <- ggplot(data = kegg.df, 
                         aes(x = Description, y = GeneRatio, fill=pvalue))+ 
        geom_bar(stat = "identity",width = 0.9)+ 
        coord_flip()+theme_bw()+ 
        scale_x_discrete(labels = function(x) str_wrap(x,width = 80))+ 
        scale_fill_gradient(low = "#E41A1C",high = "#FFF3F0")+
        labs(x = "kegg pathways",y = "GeneRatio",title = paste0("Enriched KEGG pathways of ", gene_type))+ 
        theme(axis.title = element_text(size = 13), 
              axis.text = element_text(size = 11), 
              plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
              legend.title = element_text(size = 13), 
              legend.text = element_text(size = 11), 
              plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm")) 
      ggsave(
        filename = paste0("./image/", gene_type, "_KEGG_barplot.pdf"), 
        plot = kegg_bar,
        width = 11,            
        height = 7,            
        units = "in",          
        dpi = 300              
      )
      return(kegg.df_all)
    }
  }
}


