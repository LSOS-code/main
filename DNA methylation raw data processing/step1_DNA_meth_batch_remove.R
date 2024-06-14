##### description: Batch removal and dimensionality reduction visualization for DNA methylation data

#load packages
library(Rtsne)
library(sva)
library(ChAMP)
library(ggplot2)
library(ggpubr)

##### 01Start: Dimensionality reduction for data visualization (M90/01, M180-1/02) before batch removal #####

#Load DNA methylation data of 90-day Experiment (M90),Norm_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation)
myNorm01 <- read.table("90-day spaceflight_M90/DNA methylation/Norm_beta.txt")
SampleSheet01 <- read.csv("90-day spaceflight_M90/DNA methylation/SampleSheet.csv")

#Load DNA methylation data of 180-day Experiment (M180-1),Norm_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, DNA methylation)
myNorm02 <- read.table("180-day spaceflight_M180-1/DNA methylation/Norm_beta.txt")
SampleSheet02 <- read.csv("180-day spaceflight_M180-1/DNA methylation/SampleSheet.csv")


#Merging of normalized DNA methylation data
common_prob <- intersect(rownames(myNorm01), rownames(myNorm02))
beta_mat_all <- cbind(myNorm01[common_prob,], myNorm02[common_prob,])
range(beta_mat_all)

SampleSheet_all <- rbind(SampleSheet01, SampleSheet02)
SampleSheet_all$Sample_Name <- paste0(SampleSheet_all$Sample_Name, c(rep("_90",9), rep("_180-1",9)))
SampleSheet_all$Subject <- c(paste0(SampleSheet01$Sample_Well, "_90"), paste0(SampleSheet02$Sample_Well,"_180-1"))
SampleSheet_all$Experiment <- c(rep("M90",9), rep("M180-1",9))

#Dimensionality reduction for data visualization
set.seed(123123)
tsne_out <- Rtsne(t(beta_mat_all),pca=FALSE,perplexity=2.5,theta=0.0)
tsnes=tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2") 
tsnes=as.data.frame(tsnes)
library(ggplot2)
tsnes$Subject=SampleSheet_all$Subject
colors_1 <- c("#C1548E", "#B7A637", "#37BA5F", "#90BAFF", "#F590E8", "#73B9CE")
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))+ 
      scale_color_manual(values = colors_1) +
      theme(panel.background = element_rect(fill = "white", color = "black"))+  
	  theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray")) 
tsnes$Time=SampleSheet_all$Sample_Group
tsnes$Experiment=SampleSheet_all$Experiment
colors_2 <- c("#C1548E", "#37BA5F", "#90BAFF")
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Time,shape=Experiment))+ 
      scale_color_manual(values = colors_2) +
      theme(panel.background = element_rect(fill = "white", color = "black"))+  
	  theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray")) 
meth_tsne_all <- ggpubr::ggarrange(p1, p2, nrow = 1)
write.csv(tsnes, "source data/Gene expression tsne data_Source data for Extended Data Fig. 2_a.csv", row.names=F)
##### 01End: Dimensionality reduction for data visualization before batch removal #####


##### 02Start: Removal of individual differences and batch effect by combat #####

#Removal of individual differences
myCombat01 <- champ.runCombat(beta=as.data.frame(myNorm01),pd=SampleSheet01,batchname=c("Sample_Well"))
myCombat02 <- champ.runCombat(beta=as.data.frame(myNorm02),pd=SampleSheet02,batchname=c("Sample_Well"))

#Removal of batch effect 
noBatch_beta_mat_all <- cbind(myCombat01[common_prob,], myCombat02[common_prob,])
colnames(noBatch_beta_mat_all) <- SampleSheet_all$Sample_Name
noBatch_beta_mat_all1 <- logit2(noBatch_beta_mat_all)
noBatch_beta_mat_all1 <- ComBat(dat = noBatch_beta_mat_all1, batch = SampleSheet_all$Experiment, par.prior = TRUE)
noBatch_beta_mat_all1 = ilogit2(noBatch_beta_mat_all1)
write.table(noBatch_beta_mat_all1[,1:9], file="90-day spaceflight_M90/DNA methylation/nobatch_beta.txt",quote=F)
write.table(noBatch_beta_mat_all1[,10:18], file="180-day spaceflight_M180-1/DNA methylation/nobatch_beta.txt",quote=F)
##### 02End: Removal of individual differences and batch effect by combat #####


##### 03Start: Dimensionality reduction for data visualization after batch removal #####
set.seed(123123)
tsne_out <- Rtsne(t(noBatch_beta_mat_all1),pca=FALSE,perplexity=2.5,theta=0.0)
tsnes=tsne_out$Y
colnames(tsnes) <- c("tSNE1", "tSNE2") 
tsnes=as.data.frame(tsnes)
tsnes$Subject=SampleSheet_all$Subject
colors_1 <- c("#C1548E", "#B7A637", "#37BA5F", "#90BAFF", "#F590E8", "#73B9CE")
p1 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Subject))+ 
      scale_color_manual(values = colors_1) +
      theme(panel.background = element_rect(fill = "white", color = "black"))+  
	  theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray")) 
tsnes$Time=SampleSheet_all$Sample_Group
tsnes$Experiment=SampleSheet_all$Experiment
colors_2 <- c("#C1548E", "#37BA5F", "#90BAFF")
p2 <- ggplot(tsnes, aes(x = tSNE1, y = tSNE2))+ geom_point(aes(col=Time,shape=Experiment))+ 
      scale_color_manual(values = colors_2) +
      theme(panel.background = element_rect(fill = "white", color = "black"))+  
	  theme(panel.grid.major = element_line(color = "gray"), panel.grid.minor = element_line(color = "gray")) 

write.csv(tsnes, "source data/Gene expression tsne data_Source data for Extended Data Fig. 2_b.csv", row.names=F)
noBatch_meth_tsne_all <- ggpubr::ggarrange(p1, p2, nrow = 1)

ggpubr::ggarrange(meth_tsne_all, noBatch_meth_tsne_all, ncol=1)
##### 03End: Dimensionality reduction for data visualization after batch removal #####



