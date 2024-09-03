System requirements:

softwares:
* R (v4.2.1)

R packages:
* edgeR_0.16
* foreign_0.8.86
* ggplot2_3.5.0
* ggalluvial_0.12.5
* reshape2_1.4.4
* RColorBrewer_1.1.3
* ggpubr_0.6.0

Installation guide:
#Install R software
conda install -c r r-base=4.2.1
#Install R packages in R (v4.2.1)
install.packages("edgeR")
install.packages("foreign")
install.packages("ggplot2")
install.packages("ggalluvial")
install.packages("reshape2")
install.packages("RColorBrewer")
install.packages("ggpubr")


Demo:
Load count data of 90-day experiment (M90),count_data_symbol.txt and SampleSheet.csv which are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)(saved in "90-day spaceflight_M90/RNA-seq/").  Load count data of 180-day experiment (M180-1),count_data_symbol.txt and SampleSheet.csv which are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, RNA-seq) (saved in "180-day spaceflight_M180-1/RNA-seq/").   Load count data of 180-day experiment (M180-2),count_data_symbol.txt and SampleSheet.csv which are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-2, RNA-seq) (saved in "180-day spaceflight_M180-2/RNA-seq/"). The expected results can be found in the result/ folder. The estimated time for this task is about 5 minutes.



Instructions for use:
1. Download the datasets from the provided link and  save them in a folder with the same name.
2. Install the required R packages as outlined in the installation guide.
3. Set the directory where you saved the datasets as the working directory in R.
4. Execute the R scripts in order. Step 1 performs differential expression analysis. Step 2 counts and visualizes the number of differentially expressed genes. Step 3 calculates the log fold change (logFC) of differentially expressed genes between T3 and T1 expression levels.

Note: Make sure to adjust the file paths in the script to match the location of your downloaded datasets.