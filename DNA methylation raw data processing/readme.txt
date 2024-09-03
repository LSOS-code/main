System requirements:

softwares:
* R (v4.2.1)

R packages:
* ChAMP_2.28.0
* Rtsne_0.16
* ggplot2_3.5.0
* ggpubr_0.6.0
* sva_3.46.0

Installation guide:
#Install R software
conda install -c r r-base=4.2.1
#Install R packages in R (v4.2.1)
install.packages("ChAMP")
install.packages("Rtsne")
install.packages("ggplot2")
install.packages("ggpubr")


Demo:
Download DNA methylation rawdata in  https://www.spacelifescience.cn/search/ (saved in 90-day spaceflight_M90/DNA methylation/rawdata or 180-day spaceflight_M180-1/DNA methylation/rawdata).
Load DNA methylation data of 90-day Experiment (M90),Norm_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation) (saved in "90-day spaceflight_M90/DNA methylation/").  #Load DNA methylation data of 180-day Experiment (M180-1), Norm_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, DNA methylation) (saved in "180-day spaceflight_M180-1/DNA methylation/"). The expected results can be found in the paper. The estimated time for this task is about 40 minutes.



Instructions for use:
1. Download the datasets from the provided link and  save them in a folder with the same name.
2. Install the required R packages as outlined in the installation guide.
3. Set the directory where you saved the datasets as the working directory in R.
4. Execute the R scripts in order. Step 0 shows processing of raw DNA methylation data from Illumina 850K platform. Step 1 performs batch removal and dimensionality reduction visualization for DNA methylation data.

Note: Make sure to adjust the file paths in the script to match the location of your downloaded datasets.