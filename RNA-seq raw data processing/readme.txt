System requirements:

softwares:
* FastQC software (v0.12.1)
* trim_galore software (v0.6.10)
* STAR (v2.7.10a)
* RSEM tool (v1.3.1)
* R (v4.2.1)

R packages:
* Rtsne_0.16
* ggplot2_3.5.0
* ggpubr_0.6.0
* sva_3.46.0

Installation guide:
#install fastqc
conda install -c bioconda fastqc=0.12.1
# Install Trim Galore
curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz
tar xvzf trim_galore.tar.gz
#Install STAR
wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz
tar -xzf 2.7.10a.tar.gz
cd STAR-2.7.10a/source
make STAR
#Install RSEM
conda install -c bioconda rsem=1.3.1
#Install R software
conda install -c r r-base=4.2.1
#Install R packages in R (v4.2.1)
install.packages("Rtsne")
install.packages("ggplot2")
install.packages("ggpubr")


Demo:
step0_RNA-seq_raw_data_process.sh
Download RNA-seq rawdata in  https://www.spacelifescience.cn/search/ (saved in Data/RNA-seq/rawdata/).
The expected outputs are the count data for the 90-day experiment (M90) (count_data_symbol.txt) and the count data for the 180-day experiment (M180-1) (count_data_symbol.txt). These files are available for download at https://www.spacelifescience.cn/search/ (search for '90-day spaceflight_M90, RNA-seq'  and '180-day spaceflight_M180-1, RNA-seq'). The estimated time for this task is 6-7 hours.
Download reference genome files Homo_sapiens.GRCh38.107.gtf and Homo_sapiens.GRCh38.dna.primary_assembly.fa  from http://ftp.ensembl.org/pub/release-107 (saved in Data/RNA-seq/reference_genome/).

step1_RNA-seq_batch_remove.R
Load count data of 90-day experiment (M90),count_data_symbol.txt and SampleSheet.csv which are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, RNA-seq)(saved in "90-day spaceflight_M90/RNA-seq/").  Load count data of 180-day experiment (M180-1),count_data_symbol.txt and SampleSheet.csv which are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, RNA-seq) (saved in "180-day spaceflight_M180-1/RNA-seq/").  This script enable batch removal and dimensionality reduction for enhanced data visualization of RNA-seq data. The estimated time for this task is about 10 minutes.



Instructions for use:
1. Download the datasets from the provided link and  save them in a folder with the same name.
2. Execute the script step0_RNA-seq_raw_data_process.sh in Linux: bash step0_RNA-seq_raw_data_process.sh
2. Install the required R packages as outlined in the installation guide.
3. Set the directory where you saved the datasets as the working directory in R.
4. Run step1_RNA-seq_batch_remove.R in R-4.2.1 (Linux: Rscript step1_Diff_expression_analysis.R) to perform batch removal, dimensionality reduction, and RNA-seq data visualization.

Note: Make sure to adjust the file paths in the script to match the location of your downloaded datasets.
