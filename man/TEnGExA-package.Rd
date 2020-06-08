\name{TEnGExA-package}
\alias{TEnGExA-package}
\alias{TEnGExA}
\docType{package}
\title{

\packageTitle{TEnGExA}
}
\description{

\packageDescription{TEnGExA}

}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{TEnGExA}
While working on transcripts such as CDS, alternative splicing, circular RNAs or lncRNAs, one can annotate them to assign some biological process or functions to relate with some specific pathway or network analysis. But while dealing with data from multiple tissues it is always preferable to proceed with some tissue-specific or tissue-enriched transcripts only. 
Recently tools have been developed for tissue-enrichment analysis but not only needs technical advancement to work with but more importantly limited to specific organism or predetermined list of genes or transcripts. 
R package has been developed to perform tissue-enrichment analysis of any large number of genes with any number of tissues, irrespective of any organism provided only the read count matrix.
\packageIndices{TEnGExA}

Command

result=TEnGExA(x, fpkm_flag=0, threshold=5, tissue_num=5, min_fpkm= 1)

Parameters

########## Count matrix or FPKM matrix #######

x is count matrix and matrix as below

1st column Genes Id, 2nd column length, 3rd ... nth columns are read counts

x  is fpkm values matrix as below

1st column Genes Id, 2nd column, 3rd ... nth columns are fpkm values

######## fpkm_flag 0 for count file or 1 for FPKM values ######

fpkm_flag=1

###################### Thershold is value of fpkm for expression analysis to be considered as expressed ##########

threshold=5

########  Number of tissues need to be considered for group enrichment class #

tissue_num= 5

########  FPKM value minimum threshold defualt 1  for expression analysis #########

min_fpkm= 1

######### Output ########

The output results are in matrix form with deatils in last column

}
\author{
 \packageAuthor{TEnGExA}
 
Maintainer: \packageMaintainer{TEnGExA}
}
\references{

Jain, A. and Tuteja, G. (2019) TissueEnrich: Tissue-specific gene enrichment analysis. Bioinformatics, 35(11):1966-1967. doi: 10.1093/bioinformatics/bty890.

Pei G., Dai Y., Zhao Z. and Jia P.(2019) deTS: Tissue-Specific Enrichment Analysis to decode tissue specificity. Bioinformatics, 35 : 3842-3845. doi: 10.1093/bioinformatics/btz138.

David Angeles-Albores, Raymond Y. N. Lee, Juancarlos Chan  and Paul W. Sternberg (2016) Tissue enrichment analysis for C. elegans genomics. BMC Bioinformatics 17: 366. doi:10.1186/s12859-016-1229-9

Uhlen, M. et al. (2015) Tissue-based map of the human proteome. Science, 347: 1260419-1260419.

 KEYWORDS
 
Tissue Enrich, Tissue Enhance, Gene, Tissue
}

\keyword{TEnGExA}
\seealso{deTS}
\examples{
datafile= paste(path.package("TEnGExA"),"/exdata/sample-fpkm-matrix-1.csv",sep="")
data1 = read.csv(datafile, header = TRUE)
result=TEnGExA(data1, fpkm_flag=1, threshold=5, tissue_num= 4, min_fpkm= 1)
write.csv(result, "output1.csv")

}
