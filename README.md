TEnGExA-package {TEnGExA}	

Tissue-Enrichment Analysis of any Number of Genes and Tissues

######Description ########

While working on transcripts such as CDS, alternative splicing, circular RNAs or lncRNAs, one can annotate them to assign some biological process or functions to relate with some specific pathway or network analysis. But while dealing with data from multiple tissues it is always preferable to proceed with some tissue-specific or tissue-enriched transcripts only. Recently tools have been developed for tissue-enrichment analysis but not only needs technical advancement to work with but more importantly limited to specific organism or predetermined list of genes or transcripts. R package has been developed to perform tissue-enrichment analysis of any large number of genes with any number of tissues, irrespective of any organism provided only the read count matrix.

The Webtool is also made available at URL  http://14.139.229.202/tissue_enrich/

#######Required software and packages######

 R (http://www.r-project.org/)
 
 R > 3.6.3

########Installation #######

install.packages("devtools")

library(devtools)

install_github("ubagithub/TEnGExA")


Command

result=TEnGExA(x, fpkm_flag=0, threshold=5, tissue_num=5, min_fpkm= 1)

Parameters

########## Count matrix or FPKM matrix #######

x is count matrix and matrix as below

1st column Genes Id, 2nd column length, 3rd ... nth columns are read counts

x is fpkm values matrix as below

1st column Genes Id, 2nd column, 3rd ... nth columns are fpkm values

######## fpkm_flag 0 for count file or 1 for FPKM values ######

fpkm_flag=1

###################### Thershold is value of fpkm for expression analysis to be considered as expressed ##########

threshold=5

######## Number of tissues need to be considered for group enrichment class #

tissue_num= 5

######## FPKM value minimum threshold defualt 1 for expression analysis #########

min_fpkm= 1

######### Output ########

The output results are in matrix form with deatils in last column

####### Author(s) #########

Angadi U B, Hukum C Rawal, T. K. Mondal.

Maintainer: Angadi UB : <ub.angadi@icar.gov.in>

####### References #######


####### KEYWORDS #######

Tissue Enrich, Tissue Enhance, Gene, Tissue


###### Examples  ########

datafile= paste(path.package("TEnGExA"),"/exdata/sample-fpkm-matrix-1.csv",sep="")

data1 = read.csv(datafile, header = TRUE)

result=TEnGExA(data1, fpkm_flag=1, threshold=5, tissue_num= 4, min_fpkm= 1)

write.csv(result, "output1.csv")







