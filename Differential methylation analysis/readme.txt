System requirements:

softwares:
* R (v4.2.1)

R packages:
* ChAMP_2.28.0
* ggvenn_0.1.10
* ggplot2_3.5.0
* openxlsx_4.2.5.2
* doParallel_1.0.17

Installation guide:
#Install R software
conda install -c r r-base=4.2.1
#Install R packages in R (v4.2.1)
install.packages("ChAMP")
install.packages("ggvenn")
install.packages("ggplot2")
install.packages("openxlsx")
install.packages("doParallel")


Demo:
Load DNA methylation data of 90-day Experiment (M90),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 90-day spaceflight_M90, DNA methylation) (saved in "90-day spaceflight_M90/DNA methylation/").  #Load DNA methylation data of 180-day Experiment (M180-1),nobatch_beta.txt and SampleSheet.csv are available for download at https://www.spacelifescience.cn/search/ (search for 180-day spaceflight_M180-1, DNA methylation) (saved in "180-day spaceflight_M180-1/DNA methylation/").  The expected results can be found in https://www.spacelifescience.cn/search/detail?id=24 (Results of differential methylation analysis of DNA methylation probes) and in the result/ folder . The estimated time for this task is about 20 minutes.



Instructions for use:
1. Download the datasets from the provided link and  save them in a folder with the same name.
2. Install the required R packages as outlined in the installation guide.
3. Set the directory where you saved the datasets as the working directory in R.
4. Execute the R scripts in order. Step 1 identifies differential methylation positions (DMPs). Step 2 counts the number of DMPs and presentation of the changing trend of DMPs in T2 vs T1. Step 3 identifies differential methylation regions (DMRs). Step 4 performs DMR (differentially methylated region) overlap analysis between different experiments.

Note: Make sure to adjust the file paths in the script to match the location of your downloaded datasets.
